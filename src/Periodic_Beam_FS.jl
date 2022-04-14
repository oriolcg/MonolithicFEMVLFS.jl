module Periodic_Beam_FS

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using WriteVTK
using Parameters

export run_periodic_beam_FS
export Periodic_Beam_FS_params

@with_kw struct Periodic_Beam_FS_params
  name::String = "PeriodicBeamFS"
  n::Int = 4
  dt::Real = 0.001
  tf::Real = 1.0
  order::Int = 2
  k::Int = 10
  vtk_output = false
end

function run_periodic_beam_FS(params)

  # Unpack input parameters
  @unpack name, n, dt, tf, order, k, vtk_output = params

  # Fixed parameters
  ## Geometry
  L = 2.0*π
  H = 1.0

  ## Physics
  g = 9.81
  ρ_w = 1.0e3
  ρ_b = 1.0e2
  h_b = 1.0e-2
  λ = 2*π/ k
  ω = √(g*k*tanh(k*H))
  EI_b = ρ_b*h_b*ω^2/(k^4)# + (k/kₚ)^4 - (ω/ω₀)^2
  d₀ = ρ_b*h_b/ρ_w
  Dᵨ = EI_b/ρ_w
  η₀ = 0.01
  η(x,t) = η₀*cos(k*x[1]-ω*t)
  ϕ(x,t) = η₀*ω/k * cosh(k*x[2]) / sinh(k*H) * sin(k*x[1]-ω*t)
  η(t::Real) = x -> η(x,t)
  ϕ(t::Real) = x -> ϕ(x,t)

  ## Numerics (time discretization)
  γ_t = 0.5
  β_t = 0.25
  t₀ = 0.0
  ∂uₜ_∂u = γ_t/(β_t*dt)
  ∂uₜₜ_∂u = 1/(β_t*dt^2)
  βₕ = 0.5
  αₕ = ∂uₜ_∂u/g * (1-βₕ)/βₕ

  ## Numerics (space discretization)
  h = L/n
  γ = 10.0*order*(order-1)/h

  # Define fluid domain
  println("Defining fluid domain")
  domain = (0.0, L, 0.0, H)
  partition = (2*n,n)
  𝒯_Ω = CartesianDiscreteModel(domain,partition,isperiodic=(true,false))

  # Domain size
  Lb = π
  x₀ = 0.0
  xb₀ = 0.5π
  xb₁ = xb₀ + Lb

  # Labelling
  labels_Ω = get_face_labeling(𝒯_Ω)
  add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)
  # Triangulations
  Ω = Interior(𝒯_Ω)
  Γ = Boundary(𝒯_Ω,tags="surface")

  # Auxiliar functions
  function is_beam(xs) # Check if an element is inside the beam
    n = length(xs)
    x = (1/n)*sum(xs)
    (xb₀ <= x[1] <= xb₁ ) * ( x[2] ≈ H)
  end
  function is_beam_boundary(xs) # Check if an element is on the beam boundary
    is_on_xb₀ = [x[1]≈xb₀ for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
    is_on_xb₁ = [x[1]≈xb₁ for x in xs]
    element_on_xb₀ = minimum(is_on_xb₀) # Boolean with "true" if at least one entry is true, "false" otherwise.
    element_on_xb₁ = minimum(is_on_xb₁)
    element_on_xb₀ | element_on_xb₁ # Return "true" if any of the two cases is true
  end

  # Beam triangulations
  xΓ = get_cell_coordinates(Γ)
  Γb_to_Γ_mask = lazy_map(is_beam,xΓ)
  Γb_to_Γ = findall(Γb_to_Γ_mask)
  Γf_to_Γ = findall(!,Γb_to_Γ_mask)
  Γb = Triangulation(Γ,Γb_to_Γ)
  Γfs = Triangulation(Γ,Γf_to_Γ)
  Λb = Skeleton(Γb)


  if vtk_output == true
    filename = "data/VTKOutput/5-1-4-periodic-beam-free-surface/"*name
    writevtk(Ω,filename*"_O")
    writevtk(Γ,filename*"_G")
    writevtk(Γb,filename*"_Gb")
    writevtk(Γfs,filename*"_Gfs")
    writevtk(Λb,filename*"_L")
  end

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)
  dΓb = Measure(Γb,degree)
  dΓfs = Measure(Γfs,degree)
  dΛb = Measure(Λb,degree)

  # Normals
  nΛb = get_normal_vector(Λb)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1)
  V_Γfs = TestFESpace(Γfs, reffe, conformity=:H1)
  V_Γb = TestFESpace(Γb, reffe, conformity=:H1)
  U_Ω = TransientTrialFESpace(V_Ω)
  U_Γfs = TransientTrialFESpace(V_Γfs)
  U_Γb = TransientTrialFESpace(V_Γb)
  X = TransientMultiFieldFESpace([U_Ω,U_Γfs,U_Γb])
  Y = MultiFieldFESpace([V_Ω,V_Γfs,V_Γb])

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  m((ϕₜₜ,κₜₜ,ηₜₜ),(w,u,v)) = ∫( d₀*ηₜₜ*v )dΓb
  c((ϕₜ,κₜ,ηₜ),(w,u,v)) = ∫( βₕ*ϕₜ*(u + αₕ*w) - κₜ*w )dΓfs +
                    ∫( ϕₜ*v - ηₜ*w )dΓb
  a((ϕ,κ,η),(w,u,v)) =  ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
                    ∫(  βₕ*(u + αₕ*w)*(g*κ) )dΓfs   +
                    ∫(  ( v*(g*η) + Dᵨ*Δ(v)*Δ(η) ) )dΓb +
                    ∫(  Dᵨ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb) + γ*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) )dΛb
  b((w,v)) =  ∫( 0.0 * w )dΩ
  op = TransientConstantFEOperator(m,c,a,b,X,Y)

  # Solver
  ls = LUSolver()
  ode_solver = Newmark(ls,dt,γ_t,β_t)

  # Initial solution
  x₀ = interpolate_everywhere([ϕ(0.0),η(0.0),η(0.0)],X(0.0))
  v₀ = interpolate_everywhere([∂t(ϕ)(0.0),∂t(η)(0.0),∂t(η)(0.0)],X(0.0))
  a₀ = interpolate_everywhere([∂tt(ϕ)(0.0),∂tt(η)(0.0),∂tt(η)(0.0)],X(0.0))

  # Solution
  xₜ = solve(ode_solver,op,(x₀,v₀,a₀),t₀,tf)

  # Auxiliar functions
  l2_Ω(x) = √(∑( ∫( x⋅x )dΩ ))
  l2_Γfs(x) = √(∑( ∫( x⋅x )dΓfs ))
  l2_Γb(x) = √(∑( ∫( x⋅x )dΓb ))

  t_global = Float64[]
  e_ϕ = Float64[]
  e_η = Float64[]
  E_kin_f = Float64[]
  E_pot_f = Float64[]
  E_kin_s = Float64[]
  E_ela_s = Float64[]
  ηₜ₀ = CellField(∂t(η)(0.0),Γb)
  ∇ϕ₀ = CellField(∇(ϕ(0.0)),Ω)
  Δη₀ = CellField(Δ(η(0.0)),Γb)
  η_0 = CellField(η(0.0),Γ)
  E_kin_s₀ = 0.25 * d₀ * ω^2 * η₀^2 * Lb
  E_kin_f₀ =  0.25 * g * η₀^2 * L
  E_ela_s₀ = 0.25 * Dᵨ * k^4 * η₀^2 * Lb
  E_pot_f₀ = 0.25 * g * η₀^2 * L
  if vtk_output == true
    filename = "data/VTKOutput/5-1-4-periodic-beam-free-surface/"*name
    pvd_Ω = paraview_collection(filename * "_O", append=false)
    pvd_Γ = paraview_collection(filename * "_G", append=false)
  end

  global ηₙ = x₀[2]
  global ηₙ_fv = get_free_dof_values(ηₙ)
  for ((ϕₕ,κₕ,ηₕ),tₙ) in xₜ
    push!(e_ϕ,l2_Ω(ϕ(tₙ) - ϕₕ))
    push!(e_η,l2_Γfs(η(tₙ) - κₕ)+l2_Γb(η(tₙ) - ηₕ))
    ηₜ = (ηₕ-ηₙ)/dt
    push!(E_kin_f, 0.5*∑( ∫( ∇(ϕₕ)⋅∇(ϕₕ) )dΩ ) )
    push!(E_pot_f, 0.5*g*∑( ∫( κₕ*κₕ )dΓfs ) + 0.5*g*∑( ∫( ηₕ*ηₕ )dΓb ))
    push!(E_kin_s, 0.5*d₀*∑( ∫( ηₜ*ηₜ )dΓb ) )
    push!(E_ela_s, 0.5*Dᵨ*∑( ∫( Δ(ηₕ)*Δ(ηₕ) )dΓb ) )
    push!(t_global,tₙ)

    if vtk_output == true
      pvd_Ω[tₙ] = createvtk(Ω,filename * "_O" * "_$tₙ.vtu",cellfields = ["phi" => ϕₕ],nsubcells=10)
      pvd_Γ[tₙ] = createvtk(Γ,filename * "_G" * "_$tₙ.vtu",cellfields = ["eta" => ηₕ],nsubcells=10)
    end

    ηₙ=interpolate!(ηₕ,ηₙ_fv,U_Γb(tₙ))
  end

  if vtk_output == true
    vtk_save(pvd_Ω)
    vtk_save(pvd_Γ)
  end

  return e_ϕ, e_η, E_kin_f , E_pot_f, E_kin_s, E_ela_s, E_kin_f₀, E_kin_s₀, E_pot_f₀, E_ela_s₀, t_global
end


end
