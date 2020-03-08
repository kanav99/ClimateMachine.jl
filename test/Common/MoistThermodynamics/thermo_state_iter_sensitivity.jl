using Test
using CLIMA.MoistThermodynamics
MT = MoistThermodynamics

using CLIMA.PlanetParameters
using CLIMA.RootSolvers
using LinearAlgebra
using Plots

function plot_prop(ts, dir, filename, prop)
  TS = last.(ts)
  args(i) = 1:length(getproperty(TS[i],prop)), getproperty(TS[i],prop)
  plot(args(1)..., xlabel="sa iteration", ylabel=string(prop))
  for i in 2:length(TS)
    plot!(args(i)..., xlabel="sa iteration", ylabel=string(prop))
  end
  png(joinpath(dir,filename))
end


@testset "moist thermodynamics - thermo state convergence sensitivity" begin
  FT = Float64
  ΔT_err_max = 0
  err_max = 0
  sat_adjust_call_count = 0
  e_int, ρ, q_tot, q_pt, T, p, θ_liq_ice = MT.tested_convergence_range(3, 2, 2, FT)
  RH = relative_humidity.(T, p, e_int, q_pt)
  tol = FT(1e-1)
  maxiter = 100
  domain_dim = (length(RH), length(ρ), length(T), length(maxiter));
  @show domain_dim
  TS = Array{ThermodynamicState}(undef, domain_dim);
  SOL = Array{RootSolvers.VerboseSolutionResults}(undef, domain_dim);
  results = zeros(Float64, length(maxiter));
  maxiter = maxiter isa Array ? maxiter : [maxiter,]

  mkpath("output")
  mkpath(joinpath("output","MoistThermoAnalysis"))
  dir = joinpath("output","MoistThermoAnalysis")

  # @inbounds for method in (MT.saturation_adjustment, MT.saturation_adjustment_SecantMethod)
  @inbounds for method in (MT.saturation_adjustment_SecantMethod,)
    @inbounds for (p, maxiter_) in enumerate(maxiter)
    @inbounds for (i,e_int_) in enumerate(e_int)
    @inbounds for (k,q_tot_) in enumerate(q_tot)
    @inbounds for (j,Z) in enumerate(zip(ρ, T))
      ρ_,T_ = Z
      # @show e_int_, ρ_, q_tot_, maxiter_, tol, method
      ts_eq, sol, sa_called = PhaseEquil(e_int_, ρ_, q_tot_, maxiter_, tol, method)
      ΔT_err_max = max(abs(T_ - air_temperature(ts_eq)), ΔT_err_max)
      err_max = max(abs(sol.err),err_max)
      sat_adjust_call_count += sa_called
      TS[i,j,k,p] = ts_eq
      SOL[i,j,k,p] = sol
    end
    end
    end
    end

    SOL_max_iter = SOL[:,:,:,end]
    println("------------------------------- Pass/fail rate for $(method)")
    @show length(SOL_max_iter)
    @show sum(getproperty.(SOL_max_iter, :converged))
    @show sum(getproperty.(SOL_max_iter, :converged))/length(SOL_max_iter)
    @show 1-sum(getproperty.(SOL_max_iter, :converged))/length(SOL_max_iter)
    @show ΔT_err_max
    @show err_max
    @show sat_adjust_call_count
    @show sat_adjust_call_count/length(SOL)

    # results = [count(map(x->x.converged, SOL[:,:,:,p]))/length(SOL[:,:,:,1]) for (p, maxiter_) in enumerate(maxiter)]
    # plot(maxiter, results, xlabel="maxiter", ylabel="% converged"); png(joinpath(dir,"converged_percent"))

    # results = [sum(map(x->sum(abs.(x.err)), SOL[:,:,:,p]))/length(SOL[:,:,:,1]) for (p, maxiter_) in enumerate(maxiter)]
    # plot(maxiter, results, xlabel="maxiter", ylabel="error"); png(joinpath(dir,"err_vs_maxiter"))

    # results = [sum(map(x->sum(abs.(x.iter_performed)), SOL[:,:,:,p]))/length(SOL[:,:,:,1]) for (p, maxiter_) in enumerate(maxiter)]
    # plot(maxiter, results, xlabel="maxiter", ylabel="iter_performed"); png(joinpath(dir,"iter_performed"))

    if length(SOL_max_iter)<2*10^3
      IDX_nc = findall(map(x->!x.converged, SOL[:,:,:,end]))
      IDX_co = findall(map(x->x.converged, SOL[:,:,:,end]))

      ts_nc  = zip(TS[IDX_nc],SOL[IDX_nc])
      ts_co  = zip(TS[IDX_co],SOL[IDX_co])

      plot(1:length(ts_nc), map(ts->ts.T, first.(ts_nc)), xlabel="case, no particular order", ylabel="temperature", title="non-converged"); png(joinpath(dir,"T_non_converged_$(method)"))
      plot(1:length(ts_co), map(ts->ts.T, first.(ts_co)), xlabel="case, no particular order", ylabel="temperature", title="converged"    ); png(joinpath(dir,"T_converged_$(method)"))

      if !isempty(ts_nc)
        plot_prop(ts_nc, dir, "non_converged_root_histories_$(method)", :root_history)
        plot_prop(ts_nc, dir, "non_converged_err_histories_$(method)", :err_history)
      end
      plot_prop(ts_co, dir, "converged_root_histories_$(method)", :root_history)
      plot_prop(ts_co, dir, "converged_err_histories_$(method)", :err_history)
    end

  end


end
