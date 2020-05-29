"""
    One-moment bulk microphysics scheme, which includes:

  - terminal velocity of precipitation
  - condensation and evaporation of cloud liquid water and
    deposition and sublimation of cloud ice (relaxation to equilibrium)
  - autoconversion of cloud liquid water into rain and of cloud ice into snow
  - accretion due to collisions between categories of condensed species
  - evaporation and sublimation of hydrometeors
  - melting of snow into rain
"""
module Microphysics

using ClimateMachine.MoistThermodynamics

using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, R_v, grav, T_freeze
using CLIMAParameters.Atmos.Microphysics

const APS = AbstractParameterSet
const ACloudPS = AbstractCloudParameterSet
const APrecipPS = AbstractPrecipParameterSet
const ALPS = AbstractLiquidParameterSet
const ARPS = AbstractRainParameterSet

export terminal_velocity
export conv_q_vap_to_q_liq_ice
export conv_q_liq_to_q_rai
export accretion

"""
    v0_rai(param_set, ρ)

 - `param_set` - abstract set with earth parameters
 - `ρ` air density

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
function v0_rai(param_set::APS, rain_param_set::ARPS, ρ::FT) where {FT <: Real}

    _ρ_cloud_liq::FT = ρ_cloud_liq(param_set)
    _C_drag::FT = Microphysics.C_drag(param_set)
    _grav::FT = grav(param_set)
    _r0::FT = r0(rain_param_set)

    return sqrt(FT(8 / 3) / _C_drag * (_ρ_cloud_liq / ρ - FT(1)) * _grav * _r0)
end

"""
    unpack_params(param_set, microphysics_param_set, ρ, q_)

 - `param_set` - abstract set with earth parameters
 - `microphysics_param_set` - abstract set with microphysics parameters
 - `q_` - specific humidity
 - `ρ` - air density

Utility function that unpacks microphysics parameters.
"""
@inline function unpack_params(
    param_set::APS,
    rain_param_set::ARPS,
    ρ::FT,
    q_rai::FT,
) where {FT <: Real}
    _n0::FT = n0(rain_param_set)
    _r0::FT = r0(rain_param_set)

    _m0::FT = m0(param_set, rain_param_set)
    _me::FT = me(rain_param_set)
    _a0::FT = a0(rain_param_set)
    _ae::FT = ae(rain_param_set)
    _v0::FT = v0_rai(param_set, rain_param_set, ρ)
    _ve::FT = ve(rain_param_set)

    _χm::FT = χm(rain_param_set)
    _Δm::FT = Δm(rain_param_set)
    _χa::FT = χa(rain_param_set)
    _Δa::FT = Δa(rain_param_set)
    _χv::FT = χv(rain_param_set)
    _Δv::FT = Δv(rain_param_set)

    return (
        _n0,
        _r0,
        _m0,
        _me,
        _χm,
        _Δm,
        _a0,
        _ae,
        _χa,
        _Δa,
        _v0,
        _ve,
        _χv,
        _Δv,
    )
end
"""
    lambda(q, ρ, n0, m0, me, r0, χm, Δm)

 - `q` - specific humidity of rain, ice or snow
 - `ρ` - air density
 - `n0` - size distribution parameter
 - `m0`, `me`, `χm`, `Δm`, `r0` - mass(radius) parameters

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(
    q::FT,
    ρ::FT,
    n0::FT,
    m0::FT,
    me::FT,
    r0::FT,
    χm::FT,
    Δm::FT,
) where {FT <: Real}

    λ::FT = FT(0)

    gamma_4::FT = FT(6)

    if q > FT(0)
        λ =
            (
                χm * m0 * n0 * gamma_4 / ρ / q / r0^(me + Δm)
            )^FT(1 / (me + Δm + 1))
    end
    return λ
end

"""
    terminal_velocity(param_set, precip_param_set, ρ, q_)

 - `param_set` - abstract set with earth parameters
 - `precip_param_set` - abstract set with rain or snow parameters
 - `ρ` - air density
 - `q_` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of rain drops and snow crystals.
"""
function terminal_velocity(
    param_set::APS,
    precip_param_set::APrecipPS,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, precip_param_set, ρ, q_)

        _λ::FT = lambda(q_, ρ, _n0, _m0, _me, _r0, _χm, _Δm)

        gamma_4::FT = FT(6)
        gamma_9_2::FT = FT(11.6317283965)

        fall_w =
            _χv *
            _v0 *
            (_λ * _r0)^(-_ve - _Δv) *
            gamma_9_2 / gamma_4
    end

    return fall_w
end

"""
    conv_q_vap_to_q_liq_ice(liquid_param_set, q_sat, q)

 - `liquid_param_set` - abstract set with cloud water parameters
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.
"""
function conv_q_vap_to_q_liq_ice(
    liquid_param_set::ALPS,
    q_sat::PhasePartition{FT},
    q::PhasePartition{FT},
) where {FT <: Real}

    _τ_cond_evap::FT = τ_cond_evap(liquid_param_set)

    return (q_sat.liq - q.liq) / _τ_cond_evap
end

"""
    conv_q_liq_to_q_rai(rain_param_set, q_liq)

 - `rain_param_set` - abstract set with rain microphysics parameters
 - `q_liq` - liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Kessler (1995).
"""
function conv_q_liq_to_q_rai(rain_param_set::ARPS, q_liq::FT) where {FT <: Real}

    _τ_acnv::FT = τ_acnv(rain_param_set)
    _q_liq_threshold::FT = q_liq_threshold(rain_param_set)

    return max(FT(0), q_liq - _q_liq_threshold) / _τ_acnv
end

"""
    accretion(param_set, cloud_param_set, precip_param_set, q_clo, q_pre, ρ)

 - `param_set` - abstract set with earth parameters
 - `cloud_param_set` - abstract set with cloud water or cloud ice parameters
 - `precip_param_set` - abstract set with rain or snow parameters
 - `q_clo` - cloud water or cloud ice specific humidity
 - `q_pre` - rain water or snow specific humidity
 - `ρ` - rain water or snow specific humidity

Returns the sink of cloud water (liquid or ice) due to collisions
with precipitating water (rain or snow).
"""
function accretion(
    param_set::APS,
    cloud_param_set::ACloudPS,
    precip_param_set::APrecipPS,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, precip_param_set, ρ, q_pre)

        _λ::FT = lambda(q_pre, ρ, _n0, _m0, _me, _r0, _χm, _Δm)
        _E::FT = E(cloud_param_set, precip_param_set)

        gamma_7_2::FT = FT(3.32335097044)

        accr_rate =
            q_clo * _E * _n0 * _a0 * _v0 * _χa * _χv / _λ *
            gamma_7_2 /
            (_λ * _r0)^(_ae + _ve + _Δa + _Δv)
    end
    return accr_rate
end

end #module Microphysics.jl
