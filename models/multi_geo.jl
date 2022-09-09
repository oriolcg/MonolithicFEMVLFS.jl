using Gridap
using Gridap.Io
using GridapGmsh

to_json_file(GmshDiscreteModel("multi_geo.msh"),"multi_geo.json")
to_json_file(GmshDiscreteModel("multi_geo_coarse.msh"),"multi_geo_coarse.json")
