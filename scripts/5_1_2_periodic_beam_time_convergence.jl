using DrWatson
@quickactivate "Monolithic_FEM_VLFS"

using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta

include("../src/Periodic_Beam.jl")
using .Periodic_Beam: Periodic_Beam_params, run_periodic_beam

# Extend DrWatson functions
DrWatson.allaccess(c::Periodic_Beam_params) = (:n, :dt, :tf, :order, :k)
DrWatson.default_prefix(c::Periodic_Beam_params) = c.name

# Define Execution function
function run(case::Periodic_Beam_params)
  case_name = savename(case;digits=8)
  println("-------------")
  println("Case: ",case_name)
  e_ϕ, e_η, = run_periodic_beam(case)
  e_ϕ_i = last(e_ϕ)
  e_η_i = last(e_η)
  case_name_suffix = savename(case,"jld2";digits=8)
  file = datadir("Error", case_name_suffix)
  prefix, data, suffix = DrWatson.parse_savename(case_name_suffix, parsetypes=(Int, Float64))
  push!(data, "e_ϕ_i"=>e_ϕ_i, "e_η_i"=>e_η_i)
  @tagsave(file,data)
  return data
end

# Warm-up case
k = 1
H = 1.0
g = 9.81
ω = √(g*k*tanh(k*H))
T = 2π/ω
Δt = T/1
n = 3
order = 2
path = datadir("PeriodicBeam")
case = Periodic_Beam_params(
  name="1Warm-up",
  n=n,
  dt=Δt,
  tf=T,
  k=k,
  order=order,
  vtk_output=false
)
produce_or_load(path,case,run;digits=8)

# Element size Convergence
nelem = 64
tf = 1.0
k  = 1
order = 4
path = datadir("PeriodicBeamTimeError")
e_ϕ_n = Float64[]
e_η_n = Float64[]
factor = 1.0
for i in 0:5
  Δt = tf*factor*2.0^(-i)
  case = Periodic_Beam_params(
    name="timeConvergence",
    n=nelem,
    dt=Δt,
    tf=tf,
    k=k,
    order=order
  )
  data, file = produce_or_load(path,case,run;digits=8)
  println("dt ",Δt," e ",data["e_ϕ_i"])
  push!(e_ϕ_n,data["e_ϕ_i"])
  push!(e_η_n,data["e_η_i"])
end
plot_case = Periodic_Beam_params(
  name="timeConvergence",
  dt=Δt,
  tf=tf,
  k=k,
  n=nelem,
  order=order,
)
plotName = savename(plot_case;ignores=("dt"),digits=8)

res = collect_results(path)

println("Ploting Δt-convergence")
plt1 = plot(
  fontsize=12,
  legend=:topleft,
  legendfontsize=8
)
plt2 = plot(
  fontsize=12,
  legend=:topleft,
  legendfontsize=8
)
xlabel!(plt1,"Time step size")
xlabel!(plt2,"Time step size")
ylabel!(plt1,"Error")
ylabel!(plt2,"Error")
styles = [:dash,:dashdot,:dashdotdot]
shapes = [:square,:circle,:utriangle]
res_dt  = @linq res |> where(:order .== order, :k .== k, :n .== nelem, :tf .== tf) #|> orderby(:dt)
errors_ϕ = res_dt[!,:e_ϕ_i]
errors_η = res_dt[!,:e_η_i]
Δts = res_dt[!,:dt]
println(Δts)
println(errors_η)
plot!(plt1,
  Δts,errors_ϕ,
  xaxis=:log,
  yaxis=:log,
  shape=shapes[order-1],
  color=:blue,
  style=styles[order-1],
  msize=5,
  label=L"\|\phi-\phi_h\|"
)
plot!(plt2,
  Δts,errors_η,
  xaxis=:log,
  yaxis=:log,
  shape=shapes[order-1],
  color=:red,
  style=styles[order-1],
  msize=5,
  label=L"\|\eta-\eta_h\|"#latexstring("\|\eta_h-\eta\|")
)
rate_label = latexstring("dt^{-2}")
plot!(plt1,
  Δts,0.1*Δts.^(2),
  color=:black,
  style=styles[order-1],
  label=rate_label,
  xticks=(Δts,[string(Δt) for Δt in Δts])
)
plot!(plt2,
  Δts,0.04*Δts.^(2),
  color=:black,
  style=styles[order-1],
  label=rate_label,
  xticks=(Δts,[string(Δt) for Δt in Δts])
)

savefig(plt1,plotsdir("PeriodicBeamTimeError",plotName)*"_phi.png")
savefig(plt2,plotsdir("PeriodicBeamTimeError",plotName)*"_eta.png")
