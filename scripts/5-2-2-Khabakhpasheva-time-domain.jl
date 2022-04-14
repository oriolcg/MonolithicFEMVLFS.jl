function run_5_2_2_Khavakhpasheva_time_domain()

  # Define execution function
  function run_5_2_2(case::Khabakhpasheva_time_domain_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    t,x,ηt = run_Khabakhpasheva_time_domain(case)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("5-2-2-Khabakhpasheva-time-domain", case_name_suffix)
    prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
    push!(data,"t"=>t, "x"=>x, "ηt"=>ηt)
    save(file,data)
    return data
  end

  # Warm-up case
  nx = 10
  ny = 1
  order = 2
  path = datadir("5-2-2-Khabakhpasheva-time-domain")
  case =  Khabakhpasheva_time_domain_params(
    name="2Warm-up",
    nx=nx,
    ny=ny,
    order=order
  )
  produce_or_load(path,case,run_5_2_1;digits=8)

  # Case 1: with joint
  case = Khabakhpasheva_time_domain_params(name="xi-0")
  data, file = produce_or_load(path,case,run_5_2_1)

  # Case 2: without joint
  case = Khabakhpasheva_time_domain_params(ξ=625,name="xi-625")
  data, file = produce_or_load(path,case,run_5_2_1)

  # Gather data
  res = collect_results(path)
  @show res

  # Reference data
  Khabakhpasheva_data = CSV.File(datadir("Ref_data/Khabakhpasheva","Khabakhpasheva_with_joint.csv");header=false)
  Riyansyah_data = CSV.File(datadir("Ref_data/Khabakhpasheva","Riyansyah_with_joint.csv");header=false)
  Khabakhpasheva_woj_data = CSV.File(datadir("Ref_data/Khabakhpasheva","Khabakhpasheva_without_joint.csv");header=false)
  Riyansyah_woj_data = CSV.File(datadir("Ref_data/Khabakhpasheva","Riyansyah_without_joint.csv");header=false)

  # Plot case 1
  res1 = @linq res |> where(:ξ.==0.0, :order .== 4)
  ts1 = res1[!,:t][1]
  xs1 = res1[!,:x][1]
  η_xs1 = res1[!,:ηt][1]
  ηxps = permutedims(reshape(hcat(η_xs1...),length(xs1),length(ts1)))
  η_max = [maximum(abs.(ηxps[1000:2000,it])) for it in 1:length(xs1)]
  colors = cgrad(:Blues_9,7,categorical=true,rev=true)
  plt1 = plot(xs1,η_max,
              xlims=(0,1),
              lw=2,
              label="Monolithic CG/DG (envelope)",
              palette=:rainbow)
  for (i,it) in enumerate(1000:5:1015)
    t_it = round(ts1[it],digits=3)
    plot!(plt1,xs1,η_xs1[it],
    xlims=(0,1),
    lw=1,
    label="Monolithic CG/DG t=$t_it",
    color=colors[i])
  end
  plot!(plt1,Khabakhpasheva_data.Column1,Khabakhpasheva_data.Column2,marker="o",line=false,label="Khabakhpasheva et al.")
  plot!(plt1,Riyansyah_data.Column1,Riyansyah_data.Column2,
        ls=:dash,
        lw=2,
        label="Riyansyah et al.",
        color=:green)
  xlabel!("x/L")
  ylabel!("|η|/η₀")
  savefig(plt1, plotsdir("5-2-2-Khabakhpasheva-time-domain","with_joint"))

  # Plot case 2
  res2 = @linq res |> where(:ξ.==625.0, :order .== 4)
  ts2 = res2[!,:t][1]
  xs2 = res2[!,:x][1]
  η_xs2 = res2[!,:ηt][1]
  ηxps = permutedims(reshape(hcat(η_xs2...),length(xs2),length(ts2)))
  η_max = [maximum(abs.(ηxps[1000:2000,it])) for it in 1:length(xs2)]
  plt2 = plot(xs2,η_max,
              xlims=(0,1),
              lw=2,
              label="Monolithic CG/DG (envelope)",
              palette=:rainbow)
  for (i,it) in enumerate(1000:5:1015)
    t_it = round(ts2[it],digits=3)
    plot!(plt2,xs1,η_xs2[it],
    xlims=(0,1),
    lw=1,
    label="Monolithic CG/DG t=$t_it",
    color=colors[i])
  end
  plot!(plt2,Khabakhpasheva_woj_data.Column1,Khabakhpasheva_woj_data.Column2,marker="o",line=false,label="Khabakhpasheva et al.")
  plot!(plt2,Riyansyah_woj_data.Column1,Riyansyah_woj_data.Column2,
        ls=:dash,
        lw=2,
        label="Riyansyah et al.",
        color=:green)
  xlabel!("x/L")
  ylabel!("|η|/η₀")
  savefig(plt2, plotsdir("5-2-2-Khabakhpasheva-time-domain","without_joint"))

end
