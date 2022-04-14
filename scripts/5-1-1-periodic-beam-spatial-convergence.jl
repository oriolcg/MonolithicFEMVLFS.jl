
include("../src/Periodic_Beam.jl")
using .Periodic_Beam: Periodic_Beam_params, run_periodic_beam

# Extend DrWatson functions
DrWatson.allaccess(c::Periodic_Beam_params) = (:n, :dt, :tf, :order, :k)
DrWatson.default_prefix(c::Periodic_Beam_params) = c.name

function run_5_1_1_periodic_beam_sapatial_convergence()

  # Define Execution function
  function run_5_1_1(case::Periodic_Beam_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    e_ϕ, e_η, = run_periodic_beam(case)
    e_ϕ_i = last(e_ϕ)
    e_η_i = last(e_η)
    case_name_suffix = savename(case,"jld2";digits=10)
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
  path = datadir("5-1-1-periodic-beam-spatial-convergence")
  case = Periodic_Beam_params(
    name="11Warm-up",
    n=n,
    dt=Δt,
    tf=T,
    k=k,
    order=order,
    vtk_output=false
  )
  produce_or_load(path,case,run_5_1_1;digits=8)

  # Element size Convergence
  Δt = 1.0e-6
  tf = 1.0e-4
  k  = 15
  e_ϕ_n = Float64[]
  e_η_n = Float64[]
  for order in 2:4
    for i in 1:5
      nelem = 2^(i+1)
      case = Periodic_Beam_params(
        name="spatialConvergence",
        n=nelem,
        dt=Δt,
        tf=tf,
        k=k,
        order=order
      )
      data, file = produce_or_load(path,case,run_5_1_1;digits=8)
      push!(e_ϕ_n,data["e_ϕ_i"])
      push!(e_η_n,data["e_η_i"])
    end
  end
  plot_case = Periodic_Beam_params(
    name="spatialConvergence",
    dt=Δt,
    tf=tf,
    k=k,
  )
  plotName = savename(plot_case;ignores=("n", "order"),digits=8)

  res = collect_results(path)

  println("Ploting h-convergence")
  plt1 = plot(
    fontsize=12,
    legend=:bottomleft,
    legendfontsize=8
  )
  plt2 = plot(
    fontsize=12,
    legend=:bottomleft,
    legendfontsize=8
  )
  xlabel!(plt1,"Number of elements in x-direction")
  xlabel!(plt2,"Number of elements in x-direction")
  ylabel!(plt1,"Error")
  ylabel!(plt2,"Error")
  nelems = [2^(i+1) for i in 1:5]
  styles = [:dash,:dashdot,:dashdotdot]
  shapes = [:square,:circle,:utriangle]
  factors_plt1 = [2,5,10]
  factors_plt2 = [10,30,80]
  for (iorder,order) in enumerate(2:4)
    res_order  = @linq res |> where(:order .== order, :k .== k, :dt .== Δt, :tf .== tf) |> orderby(:n)
    local_errors_ϕ = res_order[!,:e_ϕ_i]
    local_errors_η = res_order[!,:e_η_i]
    plot!(plt1,
      nelems,local_errors_ϕ,
      xaxis=:log,
      yaxis=:log,
      shape=shapes[iorder],
      color=:blue,
      style=styles[iorder],
      msize=5,
      label="r=$(order)"
    )
    plot!(plt2,
      nelems,local_errors_η,
      xaxis=:log,
      yaxis=:log,
      shape=shapes[iorder],
      color=:red,
      style=styles[iorder],
      msize=5,
      label="r=$(order+1)"
    )
  end
  for (iorder,order) in enumerate(2:4)
    rate_label = latexstring("n^{-"*"$(order+1)"*"}")
    plot!(plt1,
      nelems,factors_plt1[iorder]*nelems.^(-float(order+1)),
      color=:black,
      style=styles[iorder],
      label=rate_label,
      xticks=(nelems,[string(2*nelem) for nelem in nelems])
    )
    plot!(plt2,
      nelems,factors_plt2[iorder]*nelems.^(-float(order+1)),
      color=:black,
      style=styles[iorder],
      label=rate_label,
      xticks=(nelems,[string(2*nelem) for nelem in nelems])
    )
  end

  savefig(plt1,plotsdir("5-1-1-periodic-beam-spatial-convergence",plotName)*"_phi.png")
  savefig(plt2,plotsdir("5-1-1-periodic-beam-spatial-convergence",plotName)*"_eta.png")

end
