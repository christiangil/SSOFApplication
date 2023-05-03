## Importing packages
# using Pkg
# Pkg.activate("NEID")
# Pkg.instantiate()

using JLD2
import StellarSpectraObservationFitting as SSOF
import SSOFApplication as SSOFA

# parsing necessary inputs
gathered_jld2_fn = ARGS[1]  # something like "C:/path/to/results/ssof_rvs.jld2"
results_save_path = ARGS[2]  # something like "C:/path/to/results/". Should have subdirectories for each order
data_save_path = ARGS[3]  # something like "C:/path/to/data/". Should have subdirectories for each order

# parsing optional inputs
orders = SSOFA.order_range(ARGS[4]; first_order=7, last_order=118)  # something like "[7,118]"
results_jld2_fn = SSOF.parse_args(5, String, "results.jld2")
error_results_jld2_fn = SSOF.parse_args(6, String, "results_boot.jld2")  # usually "results_boot.jld2" or "results_curv.jld2"
data_jld2_fn = SSOF.parse_args(6, String, "data.jld2")

# making sure the inputs make sense
@assert isdir(dirname(gathered_jld2_fn)) "gathered_jld2_fn is being put in a path that doesn't exist"
@assert isdir(results_save_path) "results_save_path for input analyses does not point to a directory"
@assert isdir(data_save_path) "data_save_path for input data does not point to a directory"

error_results_jld2_jld2_full_fns = [results_save_path*"/$order/$error_results_jld2_fn" for order in orders]
results_jld2_full_fns = [results_save_path*"/$order/$results_jld2_fn" for order in orders]
data_jld2_full_fns = [data_save_path*"/$order/$data_jld2_fn" for order in orders]

@assert any(SSOFA.isjld2.(error_results_jld2_jld2_full_fns)) "No error results .jld2 files found"
@assert any(SSOFA.isjld2.(results_jld2_full_fns)) "No results .jld2 files found"
@assert any(SSOFA.isjld2.(data_jld2_full_fns)) "No data .jld2 files found"

SSOFA.retrieve(
    gathered_jld2_fn,
    error_results_jld2_jld2_full_fns,
    results_jld2_full_fns,
    data_jld2_full_fns,
    )