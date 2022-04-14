module Liu_experiment
using DrWatson
@quickactivate "Monolithic_FEM_VLFS"

using DataFrames
using DataFramesMeta
using CSV
using Plots
using LaTeXStrings

include("../src/Liu.jl")
using .Liu: Liu_params, run_Liu

# Extend DrWatson functions
DrWatson.allaccess(c::Liu_params) = (:ω,)
DrWatson.default_prefix(c::Liu_params) = c.name

# Define execution function
function run_Liu_case(case::Liu_params)
  case_name = savename(case)
  println("-------------")
  println("Case: ",case_name)
  x,η = run_Liu(case)
  case_name_suffix = savename(case,"jld2";digits=8)
  file = datadir("Liu", case_name_suffix)
  prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
  push!(data,"x"=>x, "η"=>η)
  save(file,data)
  return data
end

# Case 1: ω=0.4
path = datadir("Liu")
case = Liu_params(ω=0.4,name="omega-04")
@show case
data, file = produce_or_load(path,case,run_Liu_case)

# Case 2: ω=0.8
case = Liu_params(ω=0.8,name="omega-08")
@show case
data, file = produce_or_load(path,case,run_Liu_case)

# Gather data
res = collect_results(path)

# Reference data
Liu_data_04 = CSV.File(datadir("Liu_data","omega_04.csv");header=false)
Liu_data_08 = CSV.File(datadir("Liu_data","omega_08.csv");header=false)

# Plot case 1
res1 = @linq res |> where(:ω .== 0.4)
xs1 = res1[!,:x][1]
η_xs1 = res1[!,:η][1]
p = sortperm(xs1)
plt1 = plot(xs1[p],η_xs1[p],
            xlims=(0,1),
            ylims=(0.6,1.4),
            lw=2,
            label="Monolithic CG/DG",
            palette=:rainbow)
plot!(plt1,Liu_data_04.Column1,Liu_data_04.Column2,marker="o",line=false,label="Liu et al.")
xlabel!("x/L")
ylabel!("|η|/η₀")
savefig(plt1, plotsdir("Liu","omega_04"))

# Plot case 2
res1 = @linq res |> where(:ω .== 0.8)
xs1 = res1[!,:x][1]
η_xs1 = res1[!,:η][1]
p = sortperm(xs1)
plt1 = plot(xs1[p],η_xs1[p],
            xlims=(0,1),
            ylims=(0,1.2),
            lw=2,
            label="Monolithic CG/DG",
            palette=:rainbow)
plot!(plt1,Liu_data_08.Column1,Liu_data_08.Column2,marker="o",line=false,label="Liu et al.")
xlabel!("x/L")
ylabel!("|η|/η₀")
savefig(plt1, plotsdir("Liu","omega_08"))

end
