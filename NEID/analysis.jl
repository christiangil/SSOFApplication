## Importing packages
#using Pkg
#Pkg.activate("NEID")
# Pkg.develop(;path=".")

import StellarSpectraObservationFitting as SSOF
import SSOFApplication as SSOFA
using Statistics
using JLD2
using Plots
using DataFrames, CSV

# parsing necessary inputs
filename_order_data = ARGS[1]  # something like "C:/path/to/data/*order*/data.jld2"
neid_drp_jld2 = ARGS[2]  # something like "C:/path/to/data/neid_pipeline.jld2"
desired_order_index = SSOF.parse_args(3, Int, 81)  # used for plots and correct regularization lengthscale
base_output_path = ARGS[4]  # something like "C:/path/to/results/*order*/"

# making sure the inputs make sense
@assert SSOFA.isjld2(filename_order_data) "filename_order_data is not a .jld2 file"
@assert SSOFA.isjld2(neid_drp_jld2) "neid_drp_jld2 is not a .jld2 file"
@assert 4 <= desired_order_index <= 120 "desired_order_index is not in the right range for NEID"
@assert isdir(base_output_path) "base_output_path for output does not point to a directory"

# optional inputs
star = SSOF.parse_args(5, String, "")  # used for plot titles
n_comp_tel, n_comp_star, use_custom_n_comp =
	SSOFA.how_many_comps(SSOF.parse_args(6, String, ""))
use_custom_n_comp = SSOF.parse_args(7, Bool, false) && use_custom_n_comp
recalc = SSOF.parse_args(8, Bool, true)
log_lm = SSOF.parse_args(9, Bool, true)
dpca = SSOF.parse_args(10, Bool, false)
use_lsf = SSOF.parse_args(11, Bool, true)
use_reg = SSOF.parse_args(12, Bool, true)

interactive = length(ARGS) == 0
if !interactive; ENV["GKSwstype"] = "100" end
save_path = base_output_path * "results.jld2"

data, times_nu, airmasses = SSOFA.get_data(filename_order_data; use_lsf=use_lsf)
times_nu .-= 2400000.5

model = SSOFA.calculate_initial_model(data;
	instrument="NEID", desired_order=desired_order_index, star=star, times=times_nu,
	n_comp_tel=n_comp_tel, n_comp_star=n_comp_star, save_fn=save_path, plots_fn=base_output_path,
	recalc=recalc, use_reg=use_reg, use_custom_n_comp=use_custom_n_comp,
	dpca=dpca, log_lm=log_lm, log_λ_gp_star=1/SSOF.SOAP_gp_params.λ,
	# log_λ_gp_tel=1/110000,
	log_λ_gp_tel=1/SSOFA.neid_temporal_gp_lsf_λ(desired_order_index),
	careful_first_step=true, speed_up=false)

SSOF.no_tellurics(model) ? opt = "frozen-tel" : opt = "adam"

mws = SSOFA.create_workspace(model, data, opt)
mkpath(base_output_path*"noreg/")
df_act = SSOFA.neid_activity_indicators(neid_drp_jld2, data)  # fix this so it doesn't rely on what is in folder only
if !mws.om.metadata[:todo][:reg_improved]
	SSOFA.neid_plots(mws, airmasses, times_nu, SSOF.rvs(mws.om), zeros(length(times_nu)), base_output_path*"noreg/", neid_drp_jld2, desired_order_index;
		display_plt=interactive, df_act=df_act, title=star);
end

SSOFA.improve_regularization!(mws; save_fn=save_path, careful_first_step=true, speed_up=false)
if !mws.om.metadata[:todo][:err_estimated]; SSOFA.improve_model!(mws, airmasses, times_nu; show_plot=interactive, save_fn=save_path, iter=500, verbose=true, careful_first_step=true, speed_up=false) end
rvs, rv_errors, tel_errors, star_errors = SSOFA.estimate_σ_curvature(mws; save_fn=base_output_path * "results_curv.jld2", save_model_fn=save_path, recalc=recalc, multithread=false)
mws.om.metadata[:todo][:err_estimated] = true
rvs_b, rv_errors_b, tel_s_b, tel_errors_b, star_s_b, star_errors_b = SSOFA.estimate_σ_bootstrap(mws; save_fn=base_output_path * "results_boot.jld2", save_model_fn=save_path, recalc_mean=true, recalc=recalc, return_holders=true)

## Plots
# SSOFA.neid_plots(mws, airmasses, times_nu, rvs_b, rv_errors_b, base_output_path, neid_drp_jld2, desired_order_index;
# 	display_plt=interactive, df_act=df_act, tel_errors=tel_errors_b, star_errors=star_errors_b, title=star);
SSOFA.neid_plots(mws, airmasses, times_nu, rvs, rv_errors, base_output_path, neid_drp_jld2, desired_order_index;
	display_plt=interactive, df_act=df_act, tel_errors=tel_errors, star_errors=star_errors, title=star);