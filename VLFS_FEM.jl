module VLFS_FEM

using DrWatson
@quickactivate "VLFS_FEM"

using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta
using CSV

export run_tests

# Include source files
include("src/Periodic_Beam.jl")
include("src/Periodic_Beam_FS.jl")
include("src/Khabakhpasheva_freq_domain.jl")
include("src/Khabakhpasheva_time_domain.jl")
include("src/Liu.jl")
include("src/Yago_freq_domain.jl")

using .Periodic_Beam: Periodic_Beam_params, run_periodic_beam
using .Periodic_Beam_FS: Periodic_Beam_FS_params, run_periodic_beam_FS
using .Khabakhpasheva_freq_domain: Khabakhpasheva_freq_domain_params, run_Khabakhpasheva_freq_domain
using .Khabakhpasheva_time_domain: Khabakhpasheva_time_domain_params, run_Khabakhpasheva_time_domain
using .Liu: Liu_params, run_Liu
using .Yago_freq_domain: Yago_freq_domain_params, run_Yago_freq_domain

# Extend DrWatson functions
DrWatson.allaccess(c::Periodic_Beam_params) = (:n, :dt, :tf, :order, :k)
DrWatson.default_prefix(c::Periodic_Beam_params) = c.name
DrWatson.allaccess(c::Periodic_Beam_FS_params) = (:n, :dt, :tf, :order, :k)
DrWatson.default_prefix(c::Periodic_Beam_FS_params) = c.name
DrWatson.allaccess(c::Khabakhpasheva_freq_domain_params) = (:nx, :ny, :order, :ξ, :vtk_output)
DrWatson.default_prefix(c::Khabakhpasheva_freq_domain_params) = c.name
DrWatson.allaccess(c::Khabakhpasheva_time_domain_params) = (:nx, :ny, :order, :ξ, :vtk_output)
DrWatson.default_prefix(c::Khabakhpasheva_time_domain_params) = c.name
DrWatson.allaccess(c::Liu_params) = (:ω,)
DrWatson.default_prefix(c::Liu_params) = c.name
DrWatson.allaccess(c::Yago_freq_domain_params) = (:nx, :ny, :nz, :order, :λfactor, :dfactor)
DrWatson.default_prefix(c::Yago_freq_domain_params) = c.name

# Include script files
include("scripts/5-1-1-periodic-beam-spatial-convergence.jl")
include("scripts/5-1-2-periodic-beam-time-convergence.jl")
include("scripts/5-1-3-periodic-beam-energy.jl")
include("scripts/5-1-4-periodic-beam-free-surface-energy.jl")
include("scripts/5-2-1-Khabakhpasheva-freq-domain.jl")
include("scripts/5-2-2-Khabakhpasheva-time-domain.jl")
include("scripts/5-3-1-Liu.jl")

function run_tests(test::String)
  if test=="all"
    run_5_1_1_periodic_beam_sapatial_convergence()
    run_5_1_2_periodic_beam_time_convergence()
    run_5_1_3_periodic_beam_energy()
    run_5_1_4_periodic_beam_free_surface_energy()
    run_5_2_1_Khavakhpasheva_freq_domain()
    run_5_2_2_Khavakhpasheva_time_domain()
    run_5_3_1_Liu()
    run_5_4_1_Yago()
  elseif test == "5-1-1" || test == "5-1-1-periodic-beam-spatial-convergence"
    run_5_1_1_periodic_beam_sapatial_convergence()
  elseif test == "5-1-2" || test == "5-1-2-periodic-beam-time-convergence"
    run_5_1_2_periodic_beam_time_convergence()
  elseif test == "5-1-3" || test == "5-1-3-periodic-beam-energy"
    run_5_1_3_periodic_beam_energy()
  elseif test == "5-1-4" || test == "5-1-4-periodic-beam-free-surface-energy"
    run_5_1_4_periodic_beam_free_surface_energy()
  elseif test == "5-2-1" || test == "5-2-1-Khavakhpasheva-freq-domain"
    run_5_2_1_Khavakhpasheva_freq_domain()
  elseif test == "5-2-2" || test == "5-2-2-Khavakhpasheva-time-domain"
    run_5_2_2_Khavakhpasheva_time_domain()
  elseif test == "5-3-1" || test == "5-3-1-Liu"
    run_5_3_1_Liu()
  elseif test == "5-4-1" || test == "5-4-1-Yago"
    run_5_4_1_Yago()
  end
end

end
