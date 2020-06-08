"""
    Diagnostics

Accumulate mean fields and covariance statistics on the computational grid.
"""
module Diagnostics

export DiagnosticsGroup,
    setup_atmos_default_diagnostics,
    setup_atmos_core_diagnostics,
    setup_dump_state_and_aux_diagnostics

using CuArrays
using Dates
using FileIO
using JLD2
using KernelAbstractions
using MPI
using OrderedCollections
using Printf
using StaticArrays
import KernelAbstractions: CPU

using ..ConfigTypes
using ..DGMethods
using ..DGMethods:
    number_state_conservative,
    vars_state_conservative,
    number_state_auxiliary,
    vars_state_auxiliary,
    vars_state_gradient_flux,
    number_state_gradient_flux
using ..Mesh.Interpolation
using ..MPIStateArrays
using ..VariableTemplates
using ..Writers

using CLIMAParameters
using CLIMAParameters.Planet: planet_radius
struct EarthParameterSet <: AbstractEarthParameterSet end

Base.@kwdef mutable struct Diagnostic_Settings
    mpicomm::MPI.Comm = MPI.COMM_WORLD
    dg::Union{Nothing, DGModel} = nothing
    Q::Union{Nothing, MPIStateArray} = nothing
    starttime::Union{Nothing, String} = nothing
    param_set::Union{Nothing, AbstractEarthParameterSet} = nothing
    output_dir::Union{Nothing, String} = nothing
end
const Settings = Diagnostic_Settings()

"""
    init(mpicomm, dg, Q, starttime, output_dir, param_set)

Initialize the diagnostics collection module -- save the parameters into
`Settings`.

TODO: make `param_set` required!
"""
function init(
    mpicomm::MPI.Comm,
    dg::DGModel,
    Q::MPIStateArray,
    starttime::String,
    output_dir::String,
    param_set::AbstractEarthParameterSet = EarthParameterSet(),
)
    Settings.mpicomm = mpicomm
    Settings.dg = dg
    Settings.Q = Q
    Settings.starttime = starttime
    Settings.param_set = param_set
    Settings.output_dir = output_dir
end

include("variables.jl")
include("helpers.jl")
include("atmos_common.jl")
include("thermo.jl")
include("groups.jl")

"""
    __init()__

Module initialization function. Currently, only fills in all currently
defined diagnostic variables.
"""
function __init__()
    setup_variables()
end

end # module Diagnostics
