module Periodic_Beam

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using GridapODEs.TransientFETools
using GridapODEs.ODETools
using WriteVTK
using Parameters

export run_periodic_beam
export Periodic_Beam_params

@with_kw struct Periodic_Beam_params
  name::String = "PeriodicBeam"
  n::Int = 4
  dt::Real = 0.001
  tf::Real = 1.0
  order::Int = 2
  k::Int = 10
  vtk_output = false
end

function run_periodic_beam(params)

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

  ## Numerics (space discretization)
  h = L/n
  γ = 10.0*order*(order+1)

  # Define fluid domain
  println("Defining fluid domain")
  domain = (0.0, L, 0.0, H)
  partition = (2*n,n)
  model_Ω = CartesianDiscreteModel(domain,partition,isperiodic=(true,false))

  # Define beam domain
  println("Defining beam domain")
  labels = get_face_labeling(model_Ω)
  add_tag_from_tags!(labels,"bottom",[1,2,5])
  add_tag_from_tags!(labels,"beam",[3,4,6])

  # Triangulations
  println("Defining triangulations")
  Ω = Interior(model_Ω)
  Γ = Boundary(model_Ω,tags="beam")
  Λ = Skeleton(Γ)
  nΛ = get_normal_vector(Λ)

  # Quadratures
  println("Defining quadratures")
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)
  dΛ = Measure(Λ,degree)

  # FE spaces
  println("Defining FE spaces")
  reffe_Ω = ReferenceFE(lagrangian,Float64,order)
  reffe_Γ = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω,reffe_Ω,conformity=:H1)
  V_Γ = TestFESpace(Γ,reffe_Γ,conformity=:H1)
  U_Ω = TransientTrialFESpace(V_Ω)
  U_Γ = TransientTrialFESpace(V_Γ)
  Y = MultiFieldFESpace([V_Ω,V_Γ])
  X = TransientMultiFieldFESpace([U_Ω,U_Γ])

  # Weak form
  m((ϕₜₜ,ηₜₜ),(w,v)) = ∫( d₀*ηₜₜ*v )dΓ
  c((ϕₜ,ηₜ),(w,v)) =  ∫( ϕₜ*v - ηₜ*w )dΓ
  a((ϕ,η),(w,v)) = ∫( ∇(ϕ)⋅∇(w) )dΩ +
                   ∫( g*η*v + Dᵨ*Δ(η)*Δ(v) )dΓ +
                   ∫( Dᵨ*( - mean(Δ(η))*jump(∇(v)⋅nΛ) -
                             jump(∇(η)⋅nΛ)*mean(Δ(v)) +
                             γ/h*jump(∇(v)⋅nΛ)*jump(∇(η)⋅nΛ) ) )dΛ
  b((w,v)) =  ∫( 0.0 * w )dΩ
  op = TransientConstantFEOperator(m,c,a,b,X,Y)

  # Solver
  ls = LUSolver()
  ode_solver = Newmark(ls,dt,γ_t,β_t)

  # Initial solution
  x₀ = interpolate_everywhere([ϕ(0.0),η(0.0)],X(0.0))
  v₀ = interpolate_everywhere([∂t(ϕ)(0.0),∂t(η)(0.0)],X(0.0))
  a₀ = interpolate_everywhere([∂tt(ϕ)(0.0),∂tt(η)(0.0)],X(0.0))

  # Solution
  xₜ = solve(ode_solver,op,(x₀,v₀,a₀),t₀,tf)

  # Auxiliar functions
  l2_Ω(x) = √(∑( ∫( x⋅x )dΩ ))
  l2_Γ(x) = √(∑( ∫( x⋅x )dΓ ))

  t_global = Float64[]
  e_ϕ = Float64[]
  e_η = Float64[]
  E_kin_f = Float64[]
  E_pot_f = Float64[]
  E_kin_s = Float64[]
  E_ela_s = Float64[]
  E_kin_f₀ = 0.25 * d₀ * ω^2 * η₀^2 * L
  E_kin_s₀ = 0.25 * g * η₀^2 * L
  E_pot_f₀ = 0.25 * Dᵨ * k^4 * η₀^2 * L
  E_ela_s₀ = 0.25 * g * η₀^2 * L
  if vtk_output == true
    filename = "data/VTKOutput/Periodic_Beam/"*name
    pvd_Ω = paraview_collection(filename * "_O", append=false)
    pvd_Γ = paraview_collection(filename * "_G", append=false)
  end

  global ηₙ = x₀[2]
  global ηₙ_fv = get_free_dof_values(ηₙ)
  for ((ϕₕ,ηₕ),tₙ) in xₜ
    push!(e_ϕ,l2_Ω(ϕ(tₙ) - ϕₕ))
    push!(e_η,l2_Γ(η(tₙ) - ηₕ))
    ηₜ = (ηₕ-ηₙ)/dt
    push!(E_kin_f, 0.5*∑( ∫( ∇(ϕₕ)⋅∇(ϕₕ) )dΩ ) )
    push!(E_pot_f, 0.5*g*∑( ∫( ηₕ*ηₕ )dΓ ) )
    push!(E_kin_s, 0.5*d₀*∑( ∫( ηₜ*ηₜ )dΓ ) )
    push!(E_ela_s, 0.5*Dᵨ*∑( ∫( Δ(ηₕ)*Δ(ηₕ) )dΓ ) )
    push!(t_global,tₙ)

    if vtk_output == true
      pvd_Ω[tₙ] = createvtk(Ω,filename * "_O" * "_$tₙ.vtu",cellfields = ["phi" => ϕₕ],nsubcells=10)
      pvd_Γ[tₙ] = createvtk(Γ,filename * "_G" * "_$tₙ.vtu",cellfields = ["eta" => ηₕ],nsubcells=10)
    end

    ηₙ=interpolate!(ηₕ,ηₙ_fv,U_Γ(tₙ))
  end

  if vtk_output == true
    vtk_save(pvd_Ω)
    vtk_save(pvd_Γ)
  end

  return e_ϕ, e_η, E_kin_f , E_pot_f, E_kin_s, E_ela_s, E_kin_f₀, E_kin_s₀, E_pot_f₀, E_ela_s₀, t_global
end


end
