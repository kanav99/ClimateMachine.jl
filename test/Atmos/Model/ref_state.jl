using CLIMAParameters
using CLIMAParameters.Planet: R_d, grav, MSLP
using StaticArrays

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Geometry
using ClimateMachine.Thermodynamics
using ClimateMachine.Atmos
using ClimateMachine.ConfigTypes
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods: vars_state_auxiliary,
    number_state_auxiliary
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.VariableTemplates

using Test
using ForwardDiff

# We should provide an interface to call all physics
# kernels in some way similar to this:
function compute_ref_state(z::FT, atmos) where {FT}

    vgeo = SArray{Tuple{3,16,3}, FT}(zeros(3,16,3)) # dummy, not used
    local_geom = LocalGeometry(Val(5), vgeo, 1, 1) # dummy, not used
    st = vars_state_auxiliary(atmos, FT)
    nst = number_state_auxiliary(atmos, FT)
    arr = MArray{Tuple{nst}, FT}(undef)
    fill!(arr, 0)
    aux = Vars{st}(arr)

    # Hack: need coord in sync with incoming z, so that
    # altitude returns correct value.
    aux.coord = @SArray FT[0, 0, z]

    # Need orientation defined, so that z
    Atmos.atmos_init_aux!(
        atmos.orientation,
        atmos,
        aux,
        local_geom,
    )
    Atmos.atmos_init_aux!(
        atmos.ref_state,
        atmos,
        aux,
        local_geom,
    )
    return aux
end

@testset "Hydrostatic reference states" begin

    FT = Float64;
    m = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        moisture = DryModel(),
        init_state_conservative = x->x,
        );

    z = collect(range(FT(0), stop = FT(25e3), length = 100))

    aux_arr = compute_ref_state.(z, Ref(m))
    T = map(x->x.ref_state.T, aux_arr)
    _R_d = FT(R_d(param_set))
    _grav = FT(grav(param_set))
    _MSLP = FT(MSLP(param_set))

    function log_p_over_MSLP(_z)
        aux = compute_ref_state(_z, m)
        return log(aux.ref_state.p / _MSLP)
    end
    ∇log_p_over_MSLP =
        _z -> ForwardDiff.derivative(log_p_over_MSLP, _z)
    T_virt_rec = -_grav ./ (_R_d .* ∇log_p_over_MSLP.(z))

    # - For reference state test:
    #   - Test that `-dp/dz = \rho g`
    #   - Test that ideal gas law holds for every constructed temp:
    #   `p = \rho R_m T`, T is `aux.ref_state.T` temperature
    #   Compute temperature from total energy and make sure it
    #   equals `aux.ref_state.T`
    #   - Test that relative humidity is consistent from thermo state RH (for RH <=1)
    # - Build driver for a battery of AtmosModel test

    ρ = map(x->x.ref_state.ρ, aux_arr)
    q_tot = map(x->x.ref_state.ρq_tot, aux_arr) ./ ρ
    q_pt = PhasePartition.(q_tot)
    R_m = gas_constant_air.(Ref(param_set), q_pt)
    # Solve for `T_virt` in: T = T_virt * R_m / _R_d
    T_virt = T .* _R_d ./ R_m
    @test all(T_virt_rec .≈ T_virt)

end

