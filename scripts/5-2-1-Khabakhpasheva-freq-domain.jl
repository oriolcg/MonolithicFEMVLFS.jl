function run_5_2_1_Khavakhpasheva_freq_domain()

  # Define execution function
  function run_5_2_1(case::Khabakhpasheva_freq_domain_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    x,η = run_Khabakhpasheva_freq_domain(case)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("5-2-1-Khabakhpasheva-freq-domain", case_name_suffix)
    prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
    push!(data,"x"=>x, "η"=>η)
    save(file,data)
    return data
  end

  # Warm-up case
  nx = 10
  ny = 1
  order = 2
  path = datadir("5-2-1-Khabakhpasheva-freq-domain")
  case =  Khabakhpasheva_freq_domain_params(
    name="2Warm-up",
    nx=nx,
    ny=ny,
    order=order
  )
  produce_or_load(path,case,run_5_2_1;digits=8)

  # Case 1: with joint
  case = Khabakhpasheva_freq_domain_params(name="xi-0")
  data, file = produce_or_load(path,case,run_5_2_1)

  # Case 2: without joint
  case = Khabakhpasheva_freq_domain_params(ξ=625,name="xi-635")
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
  xs1 = res1[!,:x]
  η_xs1 = res1[!,:η]
  plt1 = plot(xs1,η_xs1,
              xlims=(0,1),
              lw=2,
              label="Monolithic CG/DG",
              palette=:rainbow)
  plot!(plt1,Khabakhpasheva_data.Column1,Khabakhpasheva_data.Column2,marker="o",line=false,label="Khabakhpasheva et al.")
  plot!(plt1,Riyansyah_data.Column1,Riyansyah_data.Column2,
        ls=:dash,
        lw=2,
        label="Riyansyah et al.")
  xlabel!("x/L")
  ylabel!("|η|/η₀")
  savefig(plt1, plotsdir("5-2-1-Khabakhpasheva-freq-domain","with_joint"))

  # Plot case 2
  res2 = @linq res |> where(:ξ.==625.0, :order .== 4)
  xs2 = res2[!,:x]
  η_xs2 = res2[!,:η]
  plt2 = plot(xs2,η_xs2,
              xlims=(0,1),
              lw=2,
              label="Monolithic CG/DG",
              palette=:rainbow)
  plot!(plt2,Khabakhpasheva_woj_data.Column1,Khabakhpasheva_woj_data.Column2,marker="o",line=false,label="Khabakhpasheva et al.")
  plot!(plt2,Riyansyah_woj_data.Column1,Riyansyah_woj_data.Column2,
        ls=:dash,
        lw=2,
        label="Riyansyah et al.")
  xlabel!("x/L")
  ylabel!("|η|/η₀")
  savefig(plt2, plotsdir("5-2-1-Khabakhpasheva-freq-domain","without_joint"))

end
