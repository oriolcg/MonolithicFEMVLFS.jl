using Gridap
using Gridap.Io
using GridapGmsh

model=GmshDiscreteModel("multi_geo.msh")
to_json_file(model,"multi_geo.json")
to_json_file(GmshDiscreteModel("multi_geo_coarse.msh"),"multi_geo_coarse.json")
writevtk(model,"multi_geo")
