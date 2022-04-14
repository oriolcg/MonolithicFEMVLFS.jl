using Gridap
using WriteVTK

# Parameters
L = 2.0*π; H = 1.0; k = 10; η₀ = 0.01; g = 9.81
ρ_w = 1.0e3; λ = 2*π/k; ω = √(g*k*tanh(k*H))
ρ_b = 1.0e2; h_b = 1.0e-2; d₀ = ρ_b*h_b/ρ_w
Dρ = ρ_b*h_b*ω^2/(k^4)/ρ_w
dt = 0.001; γ_t = 0.5; β_t = 0.25; t₀ = 0.0; tf = 1.0
order = 2; degree = 2*order; n = 30; h = L/n
γ = 1.0*order*(order+1)

# Finite element mesh
domain = (0.0,L,0.0,H); partition = (2*n,n)
model_Ω = CartesianDiscreteModel(
  domain,partition,isperiodic=(true,false))
labels = get_face_labeling(model_Ω)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"beam",[3,4,6])

# Domains
Ω = Interior(model_Ω)
Γ = Boundary(model_Ω,tags="beam")
Λ = Skeleton(Γ); nΛ = get_normal_vector(Λ)
dΩ = Measure(Ω,degree);
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω,reffe)
V_Γ = TestFESpace(Γ,reffe)
U_Ω = TransientTrialFESpace(V_Ω)
U_Γ = TransientTrialFESpace(V_Γ)
Y = MultiFieldFESpace([V_Ω,V_Γ])
X = TransientMultiFieldFESpace([U_Ω,U_Γ])

# Weak form
m((ϕₜₜ,ηₜₜ),(w,v)) = ∫( d₀*ηₜₜ*v )dΓ
c((ϕₜ,ηₜ),(w,v)) =  ∫( ϕₜ*v - ηₜ*w )dΓ
a((ϕ,η),(w,v)) =
  ∫( ∇(ϕ)⋅∇(w) )dΩ +
  ∫( g*η*v + Dρ*Δ(η)*Δ(v) )dΓ +
  ∫( Dρ*(
    - mean(Δ(η))*jump(∇(v)⋅nΛ) -
    jump(∇(η)⋅nΛ)*mean(Δ(v)) +
    γ/h*jump(∇(v)⋅nΛ)*jump(∇(η)⋅nΛ) )
   )dΛ
b((w,v)) =  ∫( 0.0 * w )dΩ
op = TransientConstantFEOperator(m,c,a,b,X,Y)

# Initial condition
η(x,t) = η₀*cos(k*x[1]-ω*t)
ϕ(x,t) =
  η₀*ω/k*cosh(k*x[2])/sinh(k*H)*sin(k*x[1]-ω*t)
η(t::Real) = x->η(x,t); ϕ(t::Real) = x->ϕ(x,t)
x₀ = interpolate_everywhere(
  [ϕ(0.0),η(0.0)],X(0.0))
v₀ = interpolate_everywhere(
  [∂t(ϕ)(0.0),∂t(η)(0.0)],X(0.0))
a₀ = interpolate_everywhere(
  [∂tt(ϕ)(0.0),∂tt(η)(0.0)],X(0.0))

# Time stepping and Paraview output
ode_solver = Newmark(LUSolver(),dt,γ_t,β_t)
xₜ = solve(ode_solver,op,(x₀,v₀,a₀),t₀,tf)
pvd_Ω = paraview_collection("Ω",append=false)
pvd_Γ = paraview_collection("Γ",append=false)
for ((ϕₕ,ηₕ),tₙ) in xₜ
  pvd_Ω[tₙ] = createvtk(
    Ω,"Ω_$tₙ.vtu",cellfields=["phi"=>ϕₕ])
  pvd_Γ[tₙ] = createvtk(
    Γ,"Γ_$tₙ.vtu",cellfields=["eta"=>ηₕ])
end
vtk_save(pvd_Ω); vtk_save(pvd_Γ)
