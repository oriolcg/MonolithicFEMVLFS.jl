function run_5_4_1_Yago()

  # Define execution function
  function run_5_4_1(case::Yago_freq_domain_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    x,η = run_Yago_freq_domain(case)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("5-4-1-Yago", case_name_suffix)
    prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
    push!(data,"x"=>x, "η"=>η)
    @tagsave(file,data)
    return data
  end

  # Warm-up case
  path = datadir("5-4-1-Yago")
  nx = 2
  ny = 2
  nz = 1
  r = 2
  case = Yago_freq_domain_params(
    name="4warm-up",
    order=r,
    nx=nx,ny=ny,nz=nz,
    λfactor=0.4,
    vtk_output=true
  )
  data, file = produce_or_load(path,case,run_5_4_1;digits=8)

  # Case 1: λfactor=0.4
  λfactor = 0.4
  nx = 32
  ny = 4
  nz = 3
  r = 4
  dfactor=4
  case = Yago_freq_domain_params(
    name="lambda-04",
    order=r,
    nx=nx,ny=ny,nz=nz,
    λfactor=λfactor,
    dfactor=dfactor,
    vtk_output=true
  )
  data, file = produce_or_load(path,case,run_5_4_1;digits=8)

  # Case 2: λfactor=0.6
  λfactor = 0.6
  case = Yago_freq_domain_params(
    name="lambda-06",
    order=r,
    nx=nx,ny=ny,nz=nz,
    λfactor=λfactor,
    dfactor=dfactor,
    vtk_output=true
  )
  data, file = produce_or_load(path,case,run_5_4_1;digits=8)

  # Case 3: λfactor=0.8
  λfactor = 0.8
  case = Yago_freq_domain_params(
    name="lambda-08",
    order=r,
    nx=nx,ny=ny,nz=nz,
    λfactor=λfactor,
    dfactor=dfactor,
    vtk_output=true
  )
  data, file = produce_or_load(path,case,run_5_4_1;digits=8)

  # Gather data
  res = collect_results(path)
  Yago_data = []
  for i in 1:10
    α = round(0.1*i,digits=2)
    push!(Yago_data, CSV.File(datadir("Ref_data/Yago","lambda_L_$α.csv");header=false))
  end
  Fu_data = []
  for case in ["L_04" "L_06" "L_08"]
    push!(Fu_data, CSV.File(datadir("Ref_data/Fu","$case.csv");header=false))
  end

  # Plot case 1
  λfactor = 0.4
  plt1 = plot(legend=:top)
  res_λ = @linq res |> where(:order .== r, :nx .== nx, :ny .== ny, :nz .== nz, :λfactor .== λfactor, :dfactor.== dfactor )
  ηs = res_λ[!,:η][1]
  xps = res_λ[!,:x][1]
  p = sortperm(xps)
  plot!(plt1,xps[p],ηs[p],label="Monolithic CG/DG",ylims=(0,1.3),lw=2,palette=:rainbow)#,marker="o",line=false)
  plot!(plt1,-2.0.*(Fu_data[1].Column1.-0.5),Fu_data[1].Column2,ls=:dash,label="Fu et al.")#color=:black,
  plot!(plt1,Yago_data[4].Column1,Yago_data[4].Column2,marker="o",line=false,label="Yago et al.")
  savefig(plt1,plotsdir("5-4-1-Yago","L_factor_04"))

  # Plot case 2
  λfactor = 0.6
  plt1 = plot(legend=:top)
  res_λ = @linq res |> where(:order .== r, :nx .== nx, :ny .== ny, :nz .== nz, :λfactor .== λfactor, :dfactor.== dfactor )
  ηs = res_λ[!,:η][1]
  xps = res_λ[!,:x][1]
  p = sortperm(xps)
  plot!(plt1,xps[p],ηs[p],label="Monolithic CG/DG",ylims=(0,1.3),lw=2,palette=:rainbow)#,color=colors[i])
  plot!(plt1,-2.0.*(Fu_data[2].Column1.-0.5),Fu_data[2].Column2,ls=:dash,label="Fu et al.")
  plot!(plt1,Yago_data[6].Column1,Yago_data[6].Column2,marker="o",line=false,label="Yago et al.")
  savefig(plt1,plotsdir("5-4-1-Yago","L_factor_06"))


  # Plot case 3
  λfactor = 0.8
  plt1 = plot(legend=:top)
  res_λ = @linq res |> where(:order .== r, :nx .== nx, :ny .== ny, :nz .== nz, :λfactor .== λfactor, :dfactor.== dfactor )
  ηs = res_λ[!,:η][1]
  xps = res_λ[!,:x][1]
  p = sortperm(xps)
  plot!(plt1,xps[p],ηs[p],label="Monolithic CG/DG",ylims=(0,1.3),lw=2,palette=:rainbow)#,color=colors[i])
  plot!(plt1,-2.0.*(Fu_data[3].Column1.-0.5),Fu_data[3].Column2,ls=:dash,label="Fu et al.")
  plot!(plt1,Yago_data[8].Column1,Yago_data[8].Column2,marker="o",line=false,label="Yago et al.")
  savefig(plt1,plotsdir("5-4-1-Yago","L_factor_08"))

end
