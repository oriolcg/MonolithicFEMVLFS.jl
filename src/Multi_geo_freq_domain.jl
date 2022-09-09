module MultiGeo_freq_domain

using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using Parameters

export run_MultiGeo_freq_domain
export MultiGeo_freq_domain_params

@with_kw struct MultiGeo_freq_domain_params
  name::String = "MultigeoFreq"
  mesh_file::String = "models/multi_geo_coarse.json"
  order::Int = 4
  dfactor::Real = 3
  vtk_output = false
end

function run_MultiGeo_freq_domain(params)

  # Unpack input parameters
  @unpack name, mesh_file, order, dfactor, vtk_output = params

  # Fixed parameters
  ## Physics
  g = 9.81
  ρ = 1025
  E = 11.9e9
  ρb = 256.25
  ν = 0.13
  D = 4.77e11
  hb = 0.1
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
  H = 3
  λ = 1.5
  k = 2π/λ
  ω = √(g*k*tanh(k*H))
  T = 2π/ω
  η₀ = 0.01
  ηᵢₙ(x) = η₀*exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*x[3]) / sinh(k*H))*exp(im*k*x[1])
  vᵢₙ(x) = (η₀*ω)*(cosh(k*x[3]) / sinh(k*H))*exp(im*k*x[1])
  vzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])

  # Define fluid model
  println("Defining fluid model")
  model_Ω = DiscreteModelFromFile(mesh_file)

  # Create masks in Ω
  labels_Ω = get_face_labeling(model_Ω)
  Γb_mask_in_Ω = get_face_mask(labels_Ω,["plate"],2)
  Γf_mask_in_Ω = get_face_mask(labels_Ω,["freeSurface"],2)

  # Triangulations
  Ω = Interior(model_Ω)
  Γ = Boundary(model_Ω,tags=["freeSurface","plate"])
  Γᵢₙ = Boundary(model_Ω,tags="inlet")
  Γb = Boundary(model_Ω,tags="plate")
  Γf = Boundary(model_Ω,tags="freeSurface")
  Λb = Skeleton(Γb)

  ## Numerics (space discretization)
  h = get_cell_measure(Λb) #0.1
  βₕ = 0.5
  αₕ = -im*ω/g * (1-βₕ)/βₕ
  γ = 1.0*order*(order+1)/h

  # Damping method 5
  μ₀ = 2.5
  Ld = 3
  Ld₀ = 3
  xd = 14
  xd₀ = 3
  μ₁(x) = μ₀*(1.0 - cos(π/2*(x[1]-xd)/Ld)) * (x[1]>xd) + μ₀*(1-cos(π/2*(Ld₀-x[1])/Ld₀)) * (x[1]<xd₀)
  μ₂(x) = μ₁(x)*k
  ηd(x) = μ₂(x)*ηᵢₙ(x)*(x[1]<xd₀)
  ∇ₙϕd(x) = μ₁(x)*vzᵢₙ(x)*(x[1]<xd₀)

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
    filename = "data/VTKOutput/5-5-1-MultiGeo/"*name
  end

  println("Computing solution")
  (ϕₕ,κₕ,ηₕ) = xₕ

  if vtk_output == true
    writevtk(Ω,filename * "_O",cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ)],nsubcells=10)
    writevtk(Γf,filename * "_Gk",cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ)],nsubcells=10)
    writevtk(Γb,filename * "_Ge",cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ)],nsubcells=10)
  end

  return nothing
end


end
