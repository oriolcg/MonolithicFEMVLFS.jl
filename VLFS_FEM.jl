module VLFS_FEM

using DrWatson
@quickactivate "VLFS_FEM"

using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta

export run

# Include source files
include("src/Periodic_Beam.jl")
using .Periodic_Beam: Periodic_Beam_params, run_periodic_beam

# Extend DrWatson functions
DrWatson.allaccess(c::Periodic_Beam_params) = (:n, :dt, :tf, :order, :k)
DrWatson.default_prefix(c::Periodic_Beam_params) = c.name

# Include script files
include("scripts/5-1-1-periodic-beam-spatial-convergence.jl")
include("scripts/5-1-2-periodic-beam-time-convergence.jl")

function run(test::String)
  if test=="all"
    run_5_1_1_periodic_beam_sapatial_convergence()
    run_5_1_2_periodic_beam_time_convergence()
  elseif test == "5-1-1" | test == "5-1-1-periodic-beam-spatial-convergence"
    run_5_1_1_periodic_beam_sapatial_convergence()
  elseif test == "5-1-2" | test == "5-1-2-periodic-beam-time-convergence"
    run_5_1_2_periodic_beam_time_convergence()
  end
end

end
