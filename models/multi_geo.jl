using Gridap
using Gridap.Io
using GridapGmsh

to_json_file(GmshDiscreteModel("multi_geo.msh"),"multi_geo.json")