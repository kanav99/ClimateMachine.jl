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

export terminal_velocity
export conv_q_vap_to_q_liq_ice
export conv_q_liq_to_q_rai
export accretion

"""
    v0_rai(ρ)

 - `ρ` air density

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
function v0_rai(ρ::FT) where {FT <: Real}

    _ρ_cloud_liq::FT = FT(1e3)
    _C_drag::FT = FT(0.55)
    _grav::FT = FT(9.81)
    _r0::FT = FT(1e-3)

    return sqrt(FT(8 / 3) / _C_drag * (_ρ_cloud_liq / ρ - FT(1)) * _grav * _r0)
end

"""
    unpack_params(ρ, q_)

 - `q_` - specific humidity
 - `ρ` - air density

Utility function that unpacks microphysics parameters.
"""
@inline function unpack_params(
    ρ::FT,
    q_rai::FT,
) where {FT <: Real}
    _n0::FT = FT(8e6*2)
    _r0::FT = FT(1e-3)
    _m0::FT = FT(4/3. * π * 1e3) * _r0^FT(3)
    _me::FT = FT(3)
    _a0::FT = FT(π) * _r0^FT(2)
    _ae::FT = FT(2)
    _v0::FT = v0_rai(ρ)
    _ve::FT = FT(0.5)

    _χm::FT = FT(1)
    _Δm::FT = FT(0)
    _χa::FT = FT(1)
    _Δa::FT = FT(0)
    _χv::FT = FT(1)
    _Δv::FT = FT(0)

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
    terminal_velocity(ρ, q_)

 - `ρ` - air density
 - `q_` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of rain drops and snow crystals.
"""
function terminal_velocity(
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(ρ, q_)

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
    conv_q_vap_to_q_liq_ice(q_sat, q)

 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.
"""
function conv_q_vap_to_q_liq_ice(
    q_sat::PhasePartition{FT},
    q::PhasePartition{FT},
) where {FT <: Real}

    _τ_cond_evap::FT = FT(10)

    return (q_sat.liq - q.liq) / _τ_cond_evap
end

"""
    conv_q_liq_to_q_rai(q_liq)

 - `q_liq` - liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Kessler (1995).
"""
function conv_q_liq_to_q_rai(q_liq::FT) where {FT <: Real}

    _τ_acnv::FT = FT(1e3)
    _q_liq_threshold::FT = FT(5e-4)

    return max(FT(0), q_liq - _q_liq_threshold) / _τ_acnv
end

"""
    accretion(q_clo, q_pre, ρ)

 - `q_clo` - cloud water or cloud ice specific humidity
 - `q_pre` - rain water or snow specific humidity
 - `ρ` - rain water or snow specific humidity

Returns the sink of cloud water (liquid or ice) due to collisions
with precipitating water (rain or snow).
"""
function accretion(
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(ρ, q_pre)

        _λ::FT = lambda(q_pre, ρ, _n0, _m0, _me, _r0, _χm, _Δm)
        _E::FT = FT(0.8)

        gamma_7_2::FT = FT(3.32335097044)

        accr_rate =
            q_clo * _E * _n0 * _a0 * _v0 * _χa * _χv / _λ *
            gamma_7_2 /
            (_λ * _r0)^(_ae + _ve + _Δa + _Δv)
    end
    return accr_rate
end

end #module Microphysics.jl
