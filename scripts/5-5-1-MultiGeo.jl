function run_5_5_1_MultiGeo()

  # Define execution function
  function run_5_5_1(case::MultiGeo_freq_domain_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    run_MultiGeo_freq_domain(case)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("5-5-1-MultiGeo", case_name_suffix)
    prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
    @tagsave(file,data)
    return data
  end

  # Warm-up case
  path = datadir("5-5-1-MultiGeo")
  r = 2
  case = MultiGeo_freq_domain_params(
    name="4warm-up",
    mesh_file="models/multi_geo_coarse.json",
    order=r,
    vtk_output=true
  )
  data, file = produce_or_load(path,case,run_5_5_1;digits=8)

  # Case 1
  r = 4
  dfactor=4
  case = MultiGeo_freq_domain_params(
    name="dfactor-4",
    mesh_file="models/multi_geo.json",
    order=r,
    dfactor=dfactor,
    vtk_output=true
  )
  data, file = produce_or_load(path,case,run_5_5_1;digits=8)

end
