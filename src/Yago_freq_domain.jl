module Yago_freq_domain

using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using Parameters

export run_Yago_freq_domain
export Yago_freq_domain_params

@with_kw struct Yago_freq_domain_params
  name::String = "YagoFreq"
  nx::Int = 32
  ny::Int = 4
  nz::Int = 4
  order::Int = 4
  λfactor::Real = 0.4
  dfactor::Real = 3
  vtk_output = false
end

function run_Yago_freq_domain(params)

  # Unpack input parameters
  @unpack name, nx, ny, nz, order, λfactor, dfactor, vtk_output = params

  # Fixed parameters
  ## Geometry
  L = 300
  B = 60
  H = 58.5
  hb = 2.0
  nLΩ = 10
  nBΩ = 18
  LΩ = nLΩ*L
  BΩ = nBΩ*B
  xb₀ = 4.5*L
  xb₁ = xb₀ + L
  yb₀ = -B/2
  yb₁ = yb₀ + B

  ## Physics
  g = 9.81
  ρ = 1025
  E = 11.9e9
  ρb = 256.25
  ν = 0.13
  D = 4.77e11
  d₀ = ρb*hb/ρ

  δ(x,y) = ==(x,y)
  C = SymFourthOrderTensorValue{3,Float64}
  μ = E/(2*(1+ν))
  λ = ν*E/(1-ν^2)
  I = hb^3/12
  Cvals = zero(Array{Float64}(undef,36))
  for i in 1:2
    for j in 1:2
      for k in 1:2
        for l in 1:2
          Cvals[data_index(C,i,j,k,l)] = I/ρ*(μ*(δ(i,k)*δ(j,l) +δ(i,l)*δ(j,k)) + λ*(δ(i,j)*δ(k,l)))
        end
      end
    end
  end
  C = SymFourthOrderTensorValue(Cvals...)

  # Wave properties
  λ = L*λfactor
  k = 2π/λ
  ω = √(g*k*tanh(k*H))
  T = 2π/ω
  η₀ = 0.01
  ηᵢₙ(x) = η₀*exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*x[3]) / sinh(k*H))*exp(im*k*x[1])
  vᵢₙ(x) = (η₀*ω)*(cosh(k*x[3]) / sinh(k*H))*exp(im*k*x[1])
  vzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])

  ## Numerics (space discretization)
  h = LΩ/(nLΩ*nx)
  βₕ = 0.5
  αₕ = -im*ω/g * (1-βₕ)/βₕ
  γ = 1.0*order*(order+1)/h

  # Damping method 5
  μ₀ = 2.5
  Ld = dfactor*L
  Ld₀ = dfactor*L
  xd = LΩ-Ld
  xd₀ = Ld₀
  μ₁(x) = μ₀*(1.0 - cos(π/2*(x[1]-xd)/Ld)) * (x[1]>xd) + μ₀*(1-cos(π/2*(Ld₀-x[1])/Ld₀)) * (x[1]<xd₀)
  μ₂(x) = μ₁(x)*k
  ηd(x) = μ₂(x)*ηᵢₙ(x)*(x[1]<xd₀)
  ∇ₙϕd(x) = μ₁(x)*vzᵢₙ(x)*(x[1]<xd₀)

  # Define fluid model
  println("Defining fluid model")
  domain = (0.0, LΩ, -BΩ/2, BΩ/2, 0.0, H)
  partition = (nLΩ*nx,nBΩ*ny,nz)
  function f_z(x)
    if x == H
        return H
    end
    i = x / (H/nz)
    return H-H/((2.5)^i)
  end
  map(x) = VectorValue(x[1], x[2], f_z(x[3]))
  model_Ω = CartesianDiscreteModel(domain,partition,map=map)

  # Add labels to Ω
  labels_Ω = get_face_labeling(model_Ω)
  add_tag_from_tags!(labels_Ω,"surface",[22])
  add_tag_from_tags!(labels_Ω,"inlet",[25])

  # Triangulations
  Ω = Interior(model_Ω)
  Γ = Boundary(model_Ω,tags="surface")
  Γᵢₙ = Boundary(model_Ω,tags="inlet")

  # Create masks in Γ
  function is_plate(x)
    is_in = ([(xb₀ <= xm[1]) * (xm[1] <= xb₁) * (yb₀ <= xm[2]) * (xm[2] <= yb₁) for xm in x])
    minimum(is_in)
  end
  xΓ = get_cell_coordinates(Γ)
  Γb_to_Γ_mask = lazy_map(is_plate,xΓ)
  Γb_to_Γ = findall(Γb_to_Γ_mask)
  Γf_to_Γ = findall(!,Γb_to_Γ_mask)
  Γb = Triangulation(Γ,Γb_to_Γ)
  Γf = Triangulation(Γ,Γf_to_Γ)

  Λb = Skeleton(Γb)

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓb = Measure(Γb,degree)
  dΓf = Measure(Γf,degree)
  dΓᵢₙ = Measure(Γᵢₙ,degree)
  dΛb = Measure(Λb,degree)

  # Normals
  nΛb = get_normal_vector(Λb)

  # FE spaces
  println("Defining FE spaces")
  reffeη = ReferenceFE(lagrangian,Float64,order)
  reffeκ = ReferenceFE(lagrangian,Float64,order)
  reffeᵩ = ReferenceFE(lagrangian,Float64,2)
  V_Ω = TestFESpace(Ω, reffeᵩ, vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γf, reffeκ, vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γb, reffeη, vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])

  # Weak form
  println("Defining weak form")
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,0.0,1.0)
  a((ϕ,κ,η),(w,u,v)) = ∫( ∇(ϕ)⋅∇(w) )dΩ +
    ∫( βₕ*(g*κ-im*ω*ϕ)*(u + αₕ*w) + im*ω*κ*w - μ₂*κ*w + μ₁*∇ₙ(ϕ)*(u + αₕ*w) )dΓf +
    ∫( ((g-ω^2*d₀)*η-im*ω*ϕ)*v + im*ω*η*w + (∇∇(v)⊙(C⊙∇∇(η))) )dΓb +
    ∫( ( - jump(∇(v))⊙(mean(C⊙∇∇(η))⋅nΛb.⁺) - (mean((C⊙∇∇(v)))⋅nΛb.⁺) ⊙jump(∇(η)) + D/ρ*γ*jump(∇(v))⊙jump(∇(η)) ) )dΛb
  b((w,u,v)) =  ∫( vᵢₙ*w )dΓᵢₙ - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓf
  op = AffineFEOperator(a,b,X,Y)

  # Solution
  println("Defining solution")
  xₕ = solve(op)

  if vtk_output == true
    filename = "data/VTKOutput/5-4-1-Yago/"*name
  end

  # Output points
  xps = []
  indices = []
  x_η = get_cell_coordinates(Γb)
  for (ie,xie) in enumerate(x_η)
    for (in,xi) in enumerate(xie)
      if -1 <= -2*(xi[1]-(xb₀+L/2))/L <=1 && xi[2]==0.0
        push!(indices,(ie,in))
        push!(xps,-2*(xi[1]-(xb₀+L/2))/L)
      end
    end
  end

  println("Computing solution")
  (ϕₕ,κₕ,ηₕ) = xₕ

  # Store values
  η_x = get_cell_dof_values(ηₕ)
  ηxps = Float64[]
  for i in 1:length(indices)
    (ie,in) = indices[i]
    push!(ηxps,abs(η_x[ie][in])/η₀)
  end

  if vtk_output == true
    writevtk(Ω,filename * "_O",cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ)],nsubcells=10)
    writevtk(Γf,filename * "_Gk",cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ)],nsubcells=10)
    writevtk(Γb,filename * "_Ge",cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ)],nsubcells=10)
  end


  return (xps,ηxps)
end


end
