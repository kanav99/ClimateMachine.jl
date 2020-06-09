#### EDMF model kernels

include(joinpath("helper_funcs", "nondimensional_exchange_functions.jl"))
include(joinpath("helper_funcs", "lamb_smooth_minimum.jl"))
include(joinpath("helper_funcs", "compute_subdomain_statistics.jl"))
include(joinpath("closures", "entr_detr.jl"))
include(joinpath("closures", "pressure.jl"))
include(joinpath("closures", "mixing_length.jl"))
# include(joinpath("closures", "micro_phys.jl"))

function vars_state_auxiliary(m::NTuple{N, Updraft}, FT) where {N}
    return Tuple{ntuple(i->vars_state_auxiliary(m[i], FT), N)...}
end

function vars_state_auxiliary(::Updraft, FT)
    @vars(buoyancy::FT,
          updraft_top::FT,
          T::FT
          )
end

function vars_state_auxiliary(::Environment, FT)
    @vars(buoyancy::FT,
          K_eddy::FT,
          cld_frac::FT,
          )
end

function vars_state_auxiliary(m::EDMF, FT)
    @vars(environment::vars_state_auxiliary(m.environment, FT),
          updraft::vars_state_auxiliary(m.updraft, FT)
          );
end

function vars_state_conservative(::Updraft, FT)
    @vars(ρa::FT,
          ρau::SVector{3, FT},
          ρae_int::FT,
          ρaq_tot::FT,
          )
end

function vars_state_conservative(::Environment, FT)
    @vars(ρatke::FT,
          ρae_int_cv::FT,
          ρaq_tot_cv::FT,
          ρae_int_q_tot_cv::FT,
          )
end

function vars_state_conservative(m::NTuple{N,Updraft}, FT) where {N}
    return Tuple{ntuple(i->vars_state_conservative(m[i], FT), N)...}
end


function vars_state_conservative(m::EDMF, FT)
    @vars(environment::vars_state_conservative(m.environment, FT),
          updraft::vars_state_conservative(m.updraft, FT)
          );
end

function vars_state_gradient(::Updraft, FT)
    @vars()
end

function vars_state_gradient(::Environment, FT)
    @vars(e_int::FT,
          q_tot::FT,
          u::SVector{3, FT},
          tke::FT,
          e_int_cv::FT,
          q_tot_cv::FT,
          e_int_q_tot_cv::FT,
          θ_ρ::FT,
          )
end


function vars_state_gradient(m::NTuple{N,Updraft}, FT) where {N}
    return Tuple{ntuple(i->vars_state_gradient(m[i], FT), N)...}
end


function vars_state_gradient(m::EDMF, FT)
    @vars(environment::vars_state_gradient(m.environment, FT),
          updraft::vars_state_gradient(m.updraft, FT)
          );
end

vars_state_gradient_flux(::Updraft, FT) = @vars()

function vars_state_gradient_flux(m::NTuple{N,Updraft}, FT) where {N}
    return Tuple{ntuple(i->vars_state_gradient_flux(m[i], FT), N)...}
end


function vars_state_gradient_flux(::Environment, FT)
    @vars(ρak∇tke::SVector{3, FT},
          ρak∇e_int_cv::SVector{3, FT},
          ρak∇q_tot_cv::SVector{3, FT},
          ρak∇e_int_q_tot_cv::SVector{3, FT},
          ∇θ_ρ::SVector{3, FT},# used in a diagnostic equation for the mixing length
          );

end

function vars_state_gradient_flux(m::EDMF, FT)
    @vars(environment::vars_state_gradient_flux(m.environment, FT),
          updraft::vars_state_gradient_flux(m.updraft, FT)
          );
end

# Specify the initial values in `aux::Vars`, which are available in
# `init_state_conservative!`. Note that
# - this method is only called at `t=0`
# - `aux.z` and `aux.T` are available here because we've specified `z` and `T`
# in `vars_state_auxiliary`
function init_state_auxiliary!(
        m::SingleStack{FT,N},
        edmf::EDMF{FT,N},
        aux::Vars,
        geom::LocalGeometry
    ) where {FT,N}
    # Aliases:
    en_a = aux.edmf.environment
    up_a = aux.edmf.updraft

    en_a.buoyancy = FT(0)
    en_a.cld_frac = FT(0)

    for i in 1:N
        up_a[i].buoyancy = FT(0)
        up_a[i].updraft_top  = FT(0)
    end
    en_a.cld_frac = FT(0)
    # Need to define K_eddy
end;

# Specify the initial values in `state::Vars`. Note that
# - this method is only called at `t=0`
# - `state.ρcT` is available here because we've specified `ρcT` in
# `vars_state_conservative`
function init_state_conservative!(
    m::SingleStack{FT,N},
    edmf::EDMF{FT,N},
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
) where {FT,N}

    # Aliases:
    gm = state
    en = state.edmf.environment
    up = state.edmf.updraft

    # gm.ρ = aux.ref_state.ρ # added at end

    # GCM setting - Initialize the grid mean profiles of prognostic variables (ρ,e_int,q_tot,u,v,w)
    z = aux.z
    gm.ρ = FT(1)
    gm.ρu = SVector(0,0,0)
    gm.ρe_int = gm.ρ * 300000
    gm.ρq_tot = gm.ρ * FT(0)

    # a_up = m.a_updraft_initial/FT(N)
    a_up = 0.1/FT(N)
    for i in 1:N
        up[i].ρa = gm.ρ * a_up
        up[i].ρau = gm.ρu * a_up
        up[i].ρae_int = gm.ρe_int * a_up
        up[i].ρaq_tot = gm.ρq_tot * a_up
    end

    # initialize environment covariance
    en.ρae_int_cv       = FT(0)
    en.ρaq_tot_cv       = FT(0)
    en.ρae_int_q_tot_cv = FT(0)

end;

# The remaining methods, defined in this section, are called at every
# time-step in the solver by the [`BalanceLaw`](@ref
# ClimateMachine.DGMethods.BalanceLaw) framework.

# Overload `update_auxiliary_state!` to call `single_stack_nodal_update_aux!`, or
# any other auxiliary methods
function update_auxiliary_state!(
    dg::DGModel,
    m::SingleStack{FT, N},
    edmf::EDMF,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
) where {FT, N}

    # Move to backend (requires `along_stack(max, var_name::Symbol)`):
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        dg.grid,
        Q,
        vars_state_conservative(m, FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        dg.grid,
        Q,
        vars_state_auxiliary(m, FT),
    )
    z = aux_vars["z"]
    for i in 1:N
        for j in 1:length(z)
            up_i_ρa = state_vars["edmf.updraft[$i].ρa"][j]
            ρinv = 1/state_vars["ρ"][j]
            if up_i_ρa*ρinv>0
                aux_vars["edmf.updraft[$i].updraft_top"][j] = z[j]
            end
        end
        # i_updraft_top = varsindex(vars_state_conservative(m, FT), :updraft_top)
        # Q[ijk, i_updraft_top, elems] = aux_vars["edmf.updraft[$i].updraft_top"][j]
    end

    nodal_update_auxiliary_state!(edmf_stack_nodal_update_aux!, dg, m, Q, t, elems)
end;

# Compute/update all auxiliary variables at each node. Note that
# - `aux.T` is available here because we've specified `T` in
# `vars_state_auxiliary`
function edmf_stack_nodal_update_aux!(
    m::SingleStack{FT,N},
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}

    en_a = aux.edmf.environment
    up_a = aux.edmf.updraft
    gm = state
    en = state.edmf.environment
    up = state.edmf.updraft

    #   -------------  Compute buoyancies of subdomains
    ρinv = 1/gm.ρ
    _grav::FT = grav(m.param_set)
    # b_upds = 0
    # a_upds = 0
    for i in 1:N
        ts = PhaseEquil(m.param_set ,up[i].ρae_int/up[i].ρa, gm.ρ, up[i].ρaq_tot/up[i].ρa)
        ρ_i = air_density(ts)
        up_a[i].buoyancy = -_grav*(ρ_i-aux.ρ0)*ρinv
    end
    # compute the buoyancy of the environment
    en_area    = 1 - sum([up[i].ρa for i in 1:N])*ρinv
    env_e_int = (gm.ρe_int - sum( [up[i].ρae_int*ρinv for i in 1:N] ))/en_area #(gm.e_int - up[i].ρae_int*ρinv)/en_area
    env_q_tot = (gm.ρq_tot - sum( [up[i].ρaq_tot*ρinv for i in 1:N] ))/en_area #(gm.q_tot - up[i].ρaq_tot*ρinv)/en_area

    ts = PhaseEquil(m.param_set ,env_e_int, gm.ρ, env_q_tot)
    env_ρ = air_density(ts)
    env_q_liq = PhasePartition(ts).liq
    env_q_ice = PhasePartition(ts).ice
    b_env = -_grav*(env_ρ - aux.ρ0)*ρinv
    b_gm = en_area*b_env + sum([up_a[i].buoyancy for i in 1:N])
    # subtract the grid mean
    for i in 1:N
        up_a[i].buoyancy -= b_gm
    end

    #   -------------  Compute upd_top
    # YAIR - check this with Charlie
    for i in 1:N
        for j in length(up[i].ρa)
            if up[i].ρa*ρinv>FT(0) # YAIR this - FT(eps) failed
                up_a[i].updraft_top = aux.z
            end
        end
    end
end;

# Since we have second-order fluxes, we must tell `ClimateMachine` to compute
# the gradient of `ρcT`. Here, we specify how `ρcT` is computed. Note that
#  - `transform.ρcT` is available here because we've specified `ρcT` in
#  `vars_state_gradient`
function compute_gradient_argument!(
    m::SingleStack{FT,N},
    edmf::EDMF{FT,N},
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}
    # Aliases:
    up_t = transform.edmf.updraft
    en_t = transform.edmf.environment
    gm = state
    up = state.edmf.updraft
    en = state.edmf.environment

    ts = thermo_state(SingleStack, state, aux)
    en_t.θ_ρ = virtual_pottemp(ts)

    ρinv       = 1/gm.ρ
    en_area    = 1 - sum([up[i].ρa for i in 1:N])*ρinv
    en_t.u     = (gm.ρu - sum([up[i].ρau for i in 1:N]))/(en_area*gm.ρ)
    en_t.e_int = (gm.ρe_int - sum([up[i].ρae_int for i in 1:N]))//(en_area*gm.ρ)
    en_t.q_tot = (gm.ρq_tot - sum([up[i].ρaq_tot for i in 1:N]))//(en_area*gm.ρ)

    en_t.tke            = en.ρatke/(en_area*gm.ρ)
    en_t.e_int_cv       = en.ρae_int_cv/(en_area*gm.ρ)
    en_t.q_tot_cv       = en.ρaq_tot_cv/(en_area*gm.ρ)
    en_t.e_int_q_tot_cv = en.ρae_int_q_tot_cv/(en_area*gm.ρ)
end;

# Specify where in `diffusive::Vars` to store the computed gradient from
# `compute_gradient_argument!`. Note that:
#  - `diffusive.α∇ρcT` is available here because we've specified `α∇ρcT` in
#  `vars_state_gradient_flux`
#  - `∇transform.ρcT` is available here because we've specified `ρcT`  in
#  `vars_state_gradient`
function compute_gradient_flux!(
    m::SingleStack{FT,N},
    edmf::EDMF{FT,N},
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT,N}
    # Aliases:
    gm_d = diffusive
    up_d = diffusive.edmf.updraft
    en_d = diffusive.edmf.environment
    gm_∇t = ∇transform
    up_∇t = ∇transform.edmf.updraft
    en_∇t = ∇transform.edmf.environment
    gm = state
    up = state.edmf.updraft
    en = state.edmf.environment

    ρinv      = 1/gm.ρ
    en_d.∇θ_ρ = en_∇t.θ_ρ

    # negative signs here as we have a '-' sign in BL form leading to + K∂ϕ/∂z on the RHS
    # compute eddy diffusivity
    εt = MArray{Tuple{N}, FT}(zeros(FT, N))
    ε = MArray{Tuple{N}, FT}(zeros(FT, N))
    δ = MArray{Tuple{N}, FT}(zeros(FT, N))
    for i in 1:N
        ε[i], δ[i], εt[i] = entr_detr(m, m.edmf.entr_detr, state, aux, t, i)
    end

    # YAIR - bug 3 
    l      = mixing_length(m, m.edmf.mix_len, state, ∇transform, aux, t, δ, εt)
    ρa_en = (gm.ρ-sum([up[j].ρa for j in 1:N]))
    K_eddy = edmf.mix_len.c_k*l*sqrt(en.ρatke/ρa_en)
    gm_d.ρak∇e_int = ρa_en*K_eddy*en_∇t.e_int
    gm_d.ρak∇q_tot = ρa_en*K_eddy*en_∇t.q_tot
    gm_d.ρak∇u     = ρa_en*K_eddy*en_∇t.u
    # second moment env cov
    en_d.ρak∇tke            = ρa_en*K_eddy*en_∇t.tke
    en_d.ρak∇e_int_cv       = ρa_en*K_eddy*en_∇t.e_int_cv
    en_d.ρak∇q_tot_cv       = ρa_en*K_eddy*en_∇t.q_tot_cv
    en_d.ρak∇e_int_q_tot_cv = ρa_en*K_eddy*en_∇t.e_int_q_tot_cv

end;

# We have no sources, nor non-diffusive fluxes.
function source!(
    m::SingleStack{FT, N},
    edmf::EDMF{FT, N},
    source::Vars,
    state::Vars,
    ∇transform::Grad,
    aux::Vars,
    t::Real,
    direction,
) where {FT, N}

    # Aliases:
    gm = state
    en = state
    up = state.edmf.updraft
    gm_s = source
    en_s = source.edmf.environment
    up_s = source.edmf.updraft

    # grid mean sources - I think that large scale subsidence in
    #            doubly periodic domains should be applied here

    # updraft sources

    # YAIR  - these need to be defined as vectors length N - check with Charlie
    εt:: SVector{N, FT}
    ε::  SVector{N, FT}
    δ::  SVector{N, FT}
    # should be conditioned on updraft_area > minval
    a_env = 1 - sum([up[i].ρa for i in 1:N])*ρinv
    w_env = (gm.ρu[3] - sum([up[i].ρau[3] for i in 1:N]))*ρinv
    e_int_env = (gm.ρe_int - sum([up[i].ρae_int for i in 1:N]))*ρinv
    q_tot_env = (gm.ρq_tot - sum([up[i].ρaq_tot for i in 1:N]))*ρinv

    ρinv  = 1/gm.ρ
    for i in 1:N
        # get environment values for e_int, q_tot , u[3]
        env_u     = (gm.u - up[i].u*up[i].ρa*ρinv)/(1-up[i].ρa*ρinv)
        env_e_int = (gm.e_int - up[i].e_int*up[i].ρa*ρinv)/(1-up[i].ρa*ρinv)
        env_q_tot = (gm.q_tot - up[i].q_tot*up[i].ρa*ρinv)/(1-up[i].ρa*ρinv)

        # first moment sources
        ε[i], δ[i], εt[i] = entr_detr(m, m.edmf.entr_detr, state, aux, t, i)
        dpdz, dpdz_tke_i = perturbation_pressure(m, m.edmf.pressure, source, state, ∇transform, aux, t, i)

        # entrainment and detrainment
        w_i = up[i].ρu[3]*ρinv
        up_s[i].ρa      += up[i].ρa * w_i * ( ε[i]                        -  δ[i])
        up_s[i].ρau     += up[i].ρa * w_i * ((ε[i]+εt[i])*up_s[i].ρau     - (δ[i]+εt[i])*env_u)
        up_s[i].ρae_int += up[i].ρa * w_i * ((ε[i]+εt[i])*up_s[i].ρae_int - (δ[i]+εt[i])*env_e_int)
        up_s[i].ρaq_tot += up[i].ρa * w_i * ((ε[i]+εt[i])*up_s[i].ρaq_tot - (δ[i]+εt[i])*env_q_tot)

        # perturbation pressure in w equation
        up_s[i].ρau[3] += up[i].ρa * dpdz
        # microphysics sources should be applied here

        ## environment second moments:

        # pressure tke source from the i'th updraft
        en.ρatke += up[i].ρa * dpdz_tke_i

        # covariances entrinament sources from the i'th updraft
        # -- if ϕ'ψ' is tke and ϕ,ψ are both w than a factor 0.5 appears in the εt and δ terms
        # Covar_Source      +=  ρaw⋅δ⋅(ϕ_up-ϕ_en)   ⋅ (ψ_up-ψ_en) + ρaw⋅εt⋅[(ϕ_up-⟨ϕ⟩)⋅(ψ_up-ψ_en) + (ϕ_up-⟨ϕ⟩)⋅(ψ_up-ψ_en)] - ρaw⋅ε⋅ϕ'ψ'
        en.ρatke            += (up[i].ρau[3] * δ[i] * (up[i].ρau[3]/up[i].ρa - w_env)*(up[i].ρau[3]/up[i].ρa - w_env)*0.5
                            +   up[i].ρau[3] * εt[i]* w_env*(up[i].ρau[3]/up[i].ρa - gm.ρu[i]*ρinv)
                            -   up[i].ρau[3] * ε[i] * en.ρatke)
        en.ρae_int_cv       += (up[i].ρau[3] * δ[i] * (up[i].ρae_int/up[i].ρa - e_int_env)*(up[i].ρae_int/up[i].ρa - e_int_env)
                            +   up[i].ρau[3] * εt[i]* e_int_env*(up[i].ρae_int/up[i].ρa - gm.ρe_int*ρinv)*2
                            -   up[i].ρau[3] * ε[i] * en.ρae_int_cv)
        en.ρaq_tot_cv       += (up[i].ρau[3] * δ[i] * (up[i].ρaq_tot/up[i].ρa - q_tot_env)*(up[i].ρaq_tot/up[i].ρa - q_tot_env)
                            +   up[i].ρau[3] * εt[i]* q_tot_env*(up[i].ρaq_tot/up[i].ρa - gm.ρq_tot*ρinv)*2
                            -   up[i].ρau[3] * ε[i] * en.ρaq_tot_cv)
        en.ρae_int_q_tot_cv += (up[i].ρau[3] * δ[i] * (up[i].ρae_int/up[i].ρa - e_int_env) *(up[i].ρaq_tot/up[i].ρa - q_tot_env)
                            +   up[i].ρau[3] * εt[i]* e_int_env*(up[i].ρaq_tot/up[i].ρa - gm.ρq_tot*ρinv)
                            +   up[i].ρau[3] * εt[i]* q_tot_env*(up[i].ρae_int/up[i].ρa - gm.ρe_int*ρinv)
                            -   up[i].ρau[3] * ε[i] * en.ρae_int_q_tot_cv)
    end
    l      = mixing_length(m, m.edmf.mix_len, state, ∇transform, aux, t, δ, εt)
    K_eddy = m.c_k*l*sqrt(en.tke)
    Shear = en_∇t.u[1].^2 + en_∇t.u[2].^2 + en_∇t.u[3].^2 # consider scalar product of two vectors

    # second moment production from mean gradients (+ sign here as we have + S in BL form)
    #                            production from mean gradient       - Dissipation
    en.ρatke            += gm.ρ*a_env*K_eddy*Shear                   - m.MixingLengthModel.c_m*sqrt(en.ρatke*ρinv/a_env)/l*en.ρatke
    en.ρae_int_cv       += gm.ρ*a_env*K_eddy*en_d.∇e_int*en_d.∇e_int - m.MixingLengthModel.c_m*sqrt(en.ρatke*ρinv/a_env)/l*en.ρae_int_cv
    en.ρaq_tot_cv       += gm.ρ*a_env*K_eddy*en_d.∇q_tot*en_d.∇q_tot - m.MixingLengthModel.c_m*sqrt(en.ρatke*ρinv/a_env)/l*en.ρaq_tot_cv
    en.ρae_int_q_tot_cv += gm.ρ*a_env*K_eddy*en_d.∇e_int*en_d.∇q_tot - m.MixingLengthModel.c_m*sqrt(en.ρatke*ρinv/a_env)/l*en.ρae_int_q_tot_cv

    # covariance microphysics sources should be applied here


end;

# in the EDMF first order (advective) fluxes exist only in the grid mean (if <w> is nonzero) and the uprdafts
function flux_first_order!(
    m::SingleStack{FT,N},
    edmf::EDMF{FT, N},
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT,N}

    # Aliases:
    gm = state
    gm_f = flux
    up = state.edmf.updraft
    up_f = flux.edmf.updraft# `vars_state_gradient_flux`

    # positive sign here as we have a '-' sign in BL form leading to - ∂ρwϕ/∂z on the RHS
    # grid mean
    ρinv = 1/gm.ρ
    gm_f.ρ = gm.ρu
    gm_f.ρe_int = gm.ρe_int*gm.ρu*ρinv
    gm_f.ρq_tot = gm.ρq_tot*gm.ρu*ρinv

    # updrafts
    for i in 1:N
        up_f[i].ρa = up[i].ρau
        u = up[i].ρau / up[i].ρa
        up_f[i].ρau = up[i].ρau * u'
        up_f[i].ρae_int = u * up[i].ρae_int
        up_f[i].ρaq_tot = u * up[i].ρaq_tot
    end

end;

# in the EDMF second order (diffusive) fluxes exist only in the grid mean and the enviroment
function flux_second_order!(
    m::SingleStack{FT,N},
    edmf::EDMF{FT,N},
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
) where {FT,N}
    # Aliases:
    gm   = state
    up   = state.edmf.updraft
    en   = state.edmf.environment
    gm_f = flux
    up_f = flux.edmf.updraft
    gm_d = diffusive
    en_d = diffusive.edmf.environment

    ρinv = FT(1)/gm.ρ
    # flux_second_order in the grid mean is the environment eddy diffusion
    en_ρa = gm.ρ-sum([up[i].ρa for i in 1:N])

    ## we are adding the massflux term here as it is part of the total flux:
    # total flux  =    diffusive_flux  +        massflux
    #   ⟨w ⃰ ϕ ⃰ ⟩  = - a_0 K_eddy⋅∂ϕ/∂z + ∑ a_i(w_i-⟨w⟩)(ϕ_i-⟨ϕ⟩)

    e_int_massflux = sum([ up[i].ρa*ρinv*(gm.ρe_int*ρinv - up[i].ρae_int/up[i].ρa)*(gm.ρu[3]*ρinv - up[i].ρau[3]/up[i].ρa) for i in 1:N])
    q_tot_massflux = sum([ up[i].ρa*ρinv*(gm.ρq_tot*ρinv - up[i].ρaq_tot/up[i].ρa)*(gm.ρu[3]*ρinv - up[i].ρau[3]/up[i].ρa) for i in 1:N])
    u_massflux     = sum([ up[i].ρa*ρinv*(gm.ρu*ρinv - up[i].ρau/up[i].ρa)*(gm.ρu[3]*ρinv - up[i].ρau[3]/up[i].ρa) for i in 1:N])

    gm_f.ρe_int  += - gm_d.k∇ρe_int + e_int_massflux
    gm_f.ρq_tot  += - gm_d.k∇ρq_tot + q_tot_massflux
    gm_f.ρu      += - gm_d.k∇ρu     + u_massflux

    # enviroment second moment turbulent flux
    en_f.ρatke            += en_d.k∇tke
    en_f.ρae_int_cv       += en_d.k∇e_int_cv
    en_f.ρaq_tot_cv       += en_d.k∇q_tot_cv
    en_f.ρae_int_q_tot_cv += en_d.k∇e_int_q_tot_cv
end;

# ### Boundary conditions

# Second-order terms in our equations, ``∇⋅(G)`` where ``G = α∇ρcT``, are
# internally reformulated to first-order unknowns.
# Boundary conditions must be specified for all unknowns, both first-order and
# second-order unknowns which have been reformulated.

# The boundary conditions for `ρcT` (first order unknown)
function boundary_state!(
    nf,
    m::SingleStack{FT,N},
    edmf::EDMF{FT,N},
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
) where {FT,N}
    gm = state⁺
    up = state⁺.edmf.updraft
    if bctype == 1 # bottom
        # placeholder to add a function for updraft surface value
        # this function should use surface covariance in the grid mean from a corresponding function
        for i in 1:N
            upd_a_surf, upd_e_int_surf, upd_q_tot_surf  = compute_updraft_surface_BC(i)
            up[i].ρau = SVector(0,0,0)
            up[i].ρa = upd_a_surf
            up[i].ρae_int = upd_e_int_surf
            up[i].ρaq_tot = upd_q_tot_surf
        end
        # can call `env_surface_covariances` with surface values

    elseif bctype == 2 # top
        # if yes not BC on upd are needed at the top (currently set to GM)
        # if not many issues might  come up with area fraction at the upd top

        for i in 1:N
            up[i].ρau = SVector(0,0,0)
            up[i].ρa = 0.0
            up[i].ρae_int = gm.ρe_int*up[i].ρa*ρinv
            up[i].ρaq_tot = gm.ρq_tot*up[i].ρa*ρinv
        end

    end
end;

# The boundary conditions for `ρcT` are specified here for second-order
# unknowns
function boundary_state!(
    nf,
    m::SingleStack{FT,N},
    edmf::EDMF{FT,N},
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
) where {FT,N}
    gm = state⁺
    up = state⁺.edmf.updraft
    gm_d = diff⁺
    up_d = diff⁺.edmf.updraft
    if bctype == 1 # bottom
        # YAIR - I need to pass the SurfaceModel into BC and into env_surface_covariances
        tke, e_int_cv ,q_tot_cv ,e_int_q_tot_cv = env_surface_covariances(m, edmf.surface, edmf, source, state)
        en_d.ρatke = gm.ρ * area_en * tke
        en_d.ρae_int_cv = gm.ρ * area_en * e_int_cv
        en_d.ρaq_tot_cv = gm.ρ * area_en * q_tot_cv
        en_d.ρae_int_q_tot_cv = gm.ρ * area_en * e_int_q_tot_cv
    elseif bctype == 2 # top
        # for now zero flux at the top
        en_d.ρatke = -n⁻ * 0.0
        en_d.ρae_int_cv = -n⁻ * 0.0
        en_d.ρaq_tot_cv = -n⁻ * 0.0
        en_d.ρae_int_q_tot_cv = -n⁻ * 0.0
    end
end;