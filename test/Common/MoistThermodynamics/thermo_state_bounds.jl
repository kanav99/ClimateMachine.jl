using Pkg

using Test
using CLIMA.MoistThermodynamics
MT = MoistThermodynamics

using CLIMA.PlanetParameters
using CLIMA.RootSolvers
using LinearAlgebra
using Plots

@testset "moist thermodynamics - bounds" begin

  FT = Float64
  e_int, ρ, q_tot, q_pt, T, p, θ_liq_ice = MT.tested_convergence_range(3, 2, 2, FT)

  domain_dim = (length(q_tot), length(ρ), length(e_int));

  TS = Array{ThermodynamicState}(undef, domain_dim)
  SOL = Array{RootSolvers.VerboseSolutionResults}(undef, domain_dim)
  maxiter = 100
  tol = 1e-1

  for (i,e_int_) in enumerate(e_int)
  for (j,ρ_) in enumerate(ρ)
  for (k,q_tot_) in enumerate(q_tot)
    ts_eq, sol, sa_called  = PhaseEquil(e_int_, ρ_, q_tot_, maxiter, tol, MT.saturation_adjustment_SecantMethod)
    TS[i,j,k] = ts_eq
    SOL[i,j,k] = sol
  end
  end
  end

# function all_phases(T, P, n_points, q)
#   M = zeros(n_points, n_points)
#   for (i, t) in enumerate(T)
#     for (j, p) in enumerate(P)
#       M[i,j] = phase(t, p, q...)
#     end
#   end
#   return M
# end
# M = all_phases(T, P, n_points, q);
# n_ice = count(x->x==1, M)
# n_liq = count(x->x==2, M)
# n_gas = count(x->x==3, M)
# contourf(T, P, (x, y)->phase(x, y, q), color=:viridis, xlabel="T", ylabel="P")
# contourf(q_tot_range, e_int_range, (x, y)->phase(x, y, q), color=:viridis, xlabel="T", ylabel="P")

  mkpath("output")
  mkpath(joinpath("output","MoistThermoAnalysis"))
  dir = joinpath("output","MoistThermoAnalysis")
  Z = Array{FT}(map(x->x.converged, SOL[1,:,:]))
  @show size(Z)
  @show size(e_int)
  @show size(ρ)

  contourf(e_int, ρ    , Z, color=:viridis, xlabel="e_int", ylabel="ρ")
  savefig(joinpath(dir,"convergedPlaneOne.png"))
  # contourf(q_tot, e_int, Array{FT}(map(x->x.converged, SOL[:,1,:])), color=:viridis, xlabel="q_tot", ylabel="e_int")
  # savefig(joinpath(dir,"converged2.png"))
  # contourf(q_tot, ρ    , Array{FT}(map(x->x.converged, SOL[:,:,1])), color=:viridis, xlabel="q_tot", ylabel="ρ")
  # savefig(joinpath(dir,"converged3.png"))

  fun(ts) = air_pressure(ts)
  fun(ts) = air_temperature(ts)
  fun(ts) = PhasePartition(ts).liq

  # contourf(e_int, ρ    , Array{FT}(map(x->fun(x), TS[1,:,:])), color=:viridis, xlabel="e_int", ylabel="ρ")
  # contourf(q_tot, e_int, map(x->fun(x), TS[:,1,:]), color=:viridis, xlabel="q_tot", ylabel="e_int")
  # contourf(q_tot, ρ    , map(x->fun(x), TS[:,:,1]), color=:viridis, xlabel="q_tot", ylabel="ρ")


end
