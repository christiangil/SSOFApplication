## Importing packages
# using Pkg
# Pkg.activate("NEID")
# Pkg.instantiate()

import NaNMath as nm
using JLD2
import StellarSpectraObservationFitting as SSOF
import SSOFApplication as SSOFA
using CSV, DataFrames
using Statistics
using Plots
using LinearAlgebra
using LaTeXStrings

# parsing necessary inputs
reduced_jld2_fn = ARGS[1]  # something like "C:/path/to/results/ssof_rvs_reduced.jld2"
fig_filename = ARGS[2]  # something like "C:/path/to/fig/ssof_neid_rv_comparison.png"
gathered_jld2_fn = ARGS[3]  # something like "C:/path/to/results/ssof_rvs.jld2"
neid_drp_jld2 = ARGS[4]  # something like "C:/path/to/data/neid_pipeline.jld2"

# parsing optional inputs
orders = SSOFA.order_range(ARGS[5]; first_order=7, last_order=118)  # something like "[7,118]"

# making sure the inputs make sense
@assert isdir(dirname(reduced_jld2_fn)) "reduced_jld2_fn is being put in a path that doesn't exist"
@assert isdir(dirname(fig_filename)) "fig_filename is being put in a path that doesn't exist"
@assert SSOFA.isjld2(gathered_jld2_fn) "gathered_jld2_fn is not a .jld2 file"
@assert SSOFA.isjld2(neid_drp_jld2) "neid_drp_jld2 is not a .jld2 file"

# loading in inputs
@load gathered_jld2_fn n_obs times_nu airmasses n_ord rvs rvs_σ constant no_tel wavelength_range wavelength_range_star all_star_s all_tel_s
times_nu .-= 2400000.5
rvs .-= mean(rvs; dims=2)

@load neid_drp_jld2 neid_time neid_rv neid_rv_σ neid_order_rv
neid_time .-= 2400000.5
neid_rv .-= mean(neid_rv)

order2ind(selected_order::Int) = SSOF.int2ind(orders, selected_order)

# http://www.csgnetwork.com/julianmodifdateconv.html
println("MJD range of observations:")
println(times_nu[1], " ", times_nu[end])

# # where is the model constant?
# star_constant = [all_star_s[i] == [] for i in eachindex(all_star_s)]
# tel_constant = [all_tel_s[i] == [] for i in eachindex(all_tel_s)]
# all(constant .== (star_constant .&& tel_constant))

# only using best orders
med_rvs_σ = vec(mean(rvs_σ; dims=2))
rvs_std = vec(std(rvs; dims=2))
best_ords = (searchsortedfirst(wavelength_range[:, 1], 4150) + 3):(searchsortedfirst(wavelength_range[:, 2], SSOF.wavenumber_to_Å(14600)) - 3)
best_ords = [i for i in best_ords if !(i in [56,57])]  # bad wavelength cal for NEID in these orders
best_ords_i = order2ind.(best_ords)
best_ords_i = best_ords_i[.!iszero.(best_ords_i)]
std_floor = 2 * nm.median(rvs_std[best_ords_i])  # round(8 * std(neid_rv))
σ_floor = 2 * nm.median(med_rvs_σ[best_ords_i])
inds_key1 = [i for i in best_ords_i if ((rvs_std[i] < std_floor) && (med_rvs_σ[i] < σ_floor))]
println("starting with $(length(best_ords_i)) orders")
println("$(length(best_ords_i) - length(inds_key1)) orders ignored for having >$std_floor m/s RMS or >$σ_floor errors or weird analysis")
println("$(length(inds_key1)) orders used in total")
rvs_red = collect(Iterators.flatten((sum(rvs[inds_key1, :] ./ (rvs_σ[inds_key1, :] .^ 2); dims=1) ./ sum(1 ./ (rvs_σ[inds_key1, :] .^ 2); dims=1))'))
rvs_red .-= median(rvs_red)
rvs_red_σ = collect(Iterators.flatten(1 ./ sqrt.(sum(1 ./ (rvs_σ[inds_key1, :] .^ 2); dims=1)')))
# plt = SSOFA.plot_rvs_red_greedy(times_nu, rvs_red, rvs_red_σ, neid_time, neid_rv, neid_rv_σ; title=star_str, inst_str="NEID");

# filtering out orders with high rms, high σ, or high χ²)
avgχ² = vec((sum(abs2, (rvs .- rvs_red') ./ rvs_σ; dims=2))) ./ n_obs
# std_floor = 2 * nm.median(rvs_std[best_ords_i])  # round(8 * std(neid_rv))
σ_floor = min(4 * nm.median(med_rvs_σ[best_ords_i]), std_floor)
filter = (rvs_std .< std_floor) .&& (med_rvs_σ .< σ_floor)
inds = [i for i in eachindex(orders) if avgχ²[i]< 4 && filter[i]]
println("starting with $(length(orders)) orders")
println("$(length(orders) - sum(filter)) orders ignored for having >$std_floor m/s RMS or >$σ_floor errors or weird analysis")
println("$(sum(filter) - length(inds)) orders ignored being uncorrelated with best order rvs")
println("$(length(inds)) orders used in total")

# using filtered orders to get more RV precision
rvs_red_greedy = collect(Iterators.flatten((sum(rvs[inds, :] ./ (rvs_σ[inds, :] .^ 2); dims=1) ./ sum(1 ./ (rvs_σ[inds, :] .^ 2); dims=1))'))
rvs_red_greedy .-= median(rvs_red_greedy)
rvs_red_greedy_σ = collect(Iterators.flatten(1 ./ sqrt.(sum(1 ./ (rvs_σ[inds, :] .^ 2); dims=1)')))
# plt = SSOFA.plot_rvs_red_greedy(times_nu, rvs_red_greedy, rvs_red_greedy_σ, neid_time, neid_rv, neid_rv_σ; title=star_str, inst_str="NEID", lower_legend=:topright)
# SSOFA.plot_rvs_red_greedy(times_nu, rvs_red_greedy, rvs_red_greedy_σ, times_nu, rvs_red, rvs_red_σ; title=star_str, inst_str="not greedy");

# making NEID vs SSOF RV plot
inst_str="Instrument";
msw=0.5;
alpha=0.7;
upper_legend=:bottomright;
lower_legend=:topright;

rvs_red_greedy .-= mean(rvs_red_greedy)
lo = minimum(append!(neid_rv .- neid_rv_σ, rvs_red_greedy .- rvs_red_greedy_σ))
hi = maximum(append!(neid_rv .+ neid_rv_σ, rvs_red_greedy .+ rvs_red_greedy_σ))
v_ssof_str = L"v^{\texttt{SSOF}}_{\star}"
v_neid_str = L"v^{\textrm{NEID}}_{\star}"

plt = SSOFA.plot_rv(; legend=upper_legend, size=(1920,1080/2))
scatter!(plt, times_nu, rvs_red_greedy; yerror=rvs_red_greedy_σ, c=SSOFA.plt_colors[2], msc=0.4*SSOFA.plt_colors[2], label=v_ssof_str * " (RMS: $(round(std(rvs_red_greedy), digits=3)), σ: $(round(mean(rvs_red_greedy_σ), digits=3)))", alpha = alpha, msw=msw)
scatter!(plt, neid_time, neid_rv; yerror=neid_rv_σ, c=SSOFA.plt_colors[1], msc=0.4*SSOFA.plt_colors[1], label=v_neid_str * " (RMS: $(round(std(neid_rv), digits=3)), σ: $(round(mean(neid_rv_σ), digits=3)))", alpha = alpha, msw=msw, ylims=(lo-0.5-(0.25*(hi-lo)),hi+0.5))

# outputs
png(fig_filename)
@save reduced_jld2_fn rvs_red rvs_red_σ rvs_red_greedy rvs_red_greedy_σ