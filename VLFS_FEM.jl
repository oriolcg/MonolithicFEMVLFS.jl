module VLFS_FEM

using DrWatson
@quickactivate "VLFS_FEM"

using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta

export run

include("scripts/5-1-1-periodic-beam-spatial-convergence.jl")

function run(test::String)
  if test=="all"
    run_5_1_1_periodic_beam_sapatial_convergence()
  elseif test == "5-1-1" | test == "5-1-1-periodic-beam-spatial-convergence"
    run_5_1_1_periodic_beam_sapatial_convergence()
  end
end

end
