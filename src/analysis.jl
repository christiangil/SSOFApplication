## Some helpful analysis functions
import StellarSpectraObservationFitting as SSOF
using JLD2
using Statistics
import StatsBase
using Base.Threads
# using ThreadsX
# using AdvancedHMC

valid_optimizers = ["adam", "l-bfgs", "frozen-tel"]


function get_data(data_fn::String; min_pix::Int=800, use_lsf::Bool=true)
	# save_path = save_path_base * star * "/$(desired_order)/"
	@load data_fn n_obs data times_nu airmasses
	if sum(all(.!(isinf.(data.var)); dims=2)) < min_pix
		@error "quitting analysis as there is not enough useful data (<$min_pix pixels used at all times)"
	end
	data.var[data.var .<= 0] .= Inf
	if !use_lsf; data = SSOF.GenericData(data) end
	return data, times_nu, airmasses
end
function calculate_initial_model(data::SSOF.Data;
	instrument::String="None", desired_order::Int=0, star::String="None", times::AbstractVector=1:size(data.flux, 2),
	n_comp_tel::Int=5, n_comp_star::Int=5, save_fn::String="", plots_fn::String="",
	recalc::Bool=false, use_custom_n_comp::Bool=false, use_reg::Bool=true, return_full_path::Bool=false, kwargs...)

	save = save_fn!=""
	if isfile(save_fn) && !recalc
		println("using saved model at $save_fn")
		@load save_fn model
	else
		oms, ℓs, aics, bics, rv_stds, rv_stds_intra, comp2ind, n_tel_cur, n_star_cur = SSOF.calculate_initial_model(data; instrument=instrument, desired_order=desired_order, star=star, times=times,
			max_n_tel=n_comp_tel, max_n_star=n_comp_star, use_all_comps=use_custom_n_comp, return_full_path=true, kwargs...)

		# plots
		if plots_fn!=""
			test_n_comp_tel = -1:n_comp_tel
			test_n_comp_star = 0:n_comp_star
			metrics = [ℓs, aics, bics, rv_stds, rv_stds_intra]
			metrics_labels = ["ℓ", "AIC", "BIC", "RV std", "Intra-night RV std"]
			metrics_fn = ["l", "aic", "bic", "rv", "rv_intra"]

			for i in eachindex(metrics)
				if any(isfinite.(metrics[i]))
					plt = component_test_plot(metrics[i], test_n_comp_tel, test_n_comp_star, ylabel=metrics_labels[i]);
					png(plt, plots_fn * metrics_fn[i] * "_choice.png")
				end
			end
		end

		if use_custom_n_comp
			model = oms[comp2ind(n_tel_cur, n_star_cur)...]
		else
			model = oms[argmin(aics)]
		end

		if !use_reg
			SSOF.rm_regularization!(model)
			model.metadata[:todo][:reg_improved] = true
		end
		if save
			@save save_fn model
		end
		if return_full_path; return oms, ℓs, aics, bics, rv_stds, rv_stds_intra, comp2ind, n_tel_cur, n_star_cur end
	end
	return model
end

function create_workspace(model, data, opt::String)
	@assert opt in valid_optimizers
	if opt == "l-bfgs"
		mws = SSOF.OptimTelStarWorkspace(model, data)
	elseif opt == "frozen-tel"
		mws = SSOF.FrozenTelWorkspace(model, data)
	else
		mws = SSOF.TotalWorkspace(model, data)
	end
	return mws
end

function improve_regularization!(mws::SSOF.ModelWorkspace; redo::Bool=false, verbose::Bool=true, testing_ratio::Real=0.33, save_fn::String="", careful_first_step::Bool=true, speed_up::Bool=false, kwargs...)

	save = save_fn!=""

	model = mws.om
	if redo || !model.metadata[:todo][:reg_improved]  # 27 mins
		@assert 0 < testing_ratio < 1
		n_obs = size(mws.d.flux, 2)

	    SSOF.train_OrderModel!(mws; verbose=verbose, ignore_regularization=true, careful_first_step=careful_first_step, speed_up=speed_up)  # 45s
	    n_obs_test = Int(round(testing_ratio * n_obs))
	    test_start_ind = max(1, Int(round(rand() * (n_obs - n_obs_test))))
	    testing_inds = test_start_ind:test_start_ind+n_obs_test-1
	    SSOF.fit_regularization!(mws, testing_inds; kwargs...)
	    model.metadata[:todo][:reg_improved] = true
	    if save; @save save_fn model end
	end
end

function improve_model!(mws::SSOF.ModelWorkspace; verbose::Bool=true, show_plot::Bool=false, save_fn::String="", kwargs...)
	save = save_fn!=""
	model = mws.om
	SSOF.train_OrderModel!(mws; verbose=verbose, kwargs...)  # 120s
	results = SSOF.finalize_scores!(mws; verbose=verbose, kwargs...)
    if show_plot; status_plot(mws) end
    if save; @save save_fn model end
	return results
end
function improve_model!(mws::SSOF.ModelWorkspace, airmasses::AbstractVector, times::AbstractVector; show_plot::Bool=false, kwargs...)
	results = improve_model!(mws; show_plot=show_plot, kwargs...)
	if show_plot; plot_model(mws, airmasses, times) end
	return results
end


function estimate_σ_curvature(mws::SSOF.ModelWorkspace; recalc::Bool=false, save_fn::String="", save_model_fn::String="", kwargs...)
	save = save_fn!=""
	save_model = save_model_fn==""
	model = mws.om
	if recalc || !model.metadata[:todo][:err_estimated] #|| !isfile(save_fn)
	    mws.d.var[mws.d.var.==Inf] .= 0
	    data_noise = sqrt.(mws.d.var)
	    mws.d.var[mws.d.var.==0] .= Inf

		time_var_tel = SSOF.is_time_variable(model.tel)
		time_var_star = SSOF.is_time_variable(model.star)

		typeof(model) <: SSOF.OrderModelDPCA ? rvs = copy(model.rv.lm.s) : rvs = copy(model.rv)
		ℓ_rv(x) = SSOF._loss(mws.o, model, mws.d; rv=x) / 2  # factor of 2 makes curvature estimates correct (χ² -> data fit part of multivariate Gaussian)
		rvs_σ = estimate_σ_curvature_helper(rvs, ℓ_rv; param_str="rv", kwargs...)
		if typeof(model) <: SSOF.OrderModelDPCA
			rvs = vec(rvs)
			rvs .*= -SSOF.light_speed_nu
			rvs_σ .*= SSOF.light_speed_nu
		end

		if time_var_tel
			ℓ_tel(x) = (SSOF._loss(mws.o, model, mws.d; tel=vec(model.tel.lm)) + SSOF.model_s_prior(model.tel.lm.s, model.reg_tel)) / 2  # factor of 2 makes curvature estimates correct (χ² -> data fit part of multivariate Gaussian)
			tel_s_σ = estimate_σ_curvature_helper(model.tel.lm.s, ℓ_tel; param_str="tel_s", kwargs...)
		else
			tel_s_σ = nothing
		end

		if time_var_star
			ℓ_star(x) = (SSOF._loss(mws.o, model, mws.d; star=vec(model.star.lm)) + SSOF.model_s_prior(model.star.lm.s, model.reg_star)) / 2  # factor of 2 makes curvature estimates correct (χ² -> data fit part of multivariate Gaussian)
			star_s_σ = estimate_σ_curvature_helper(model.star.lm.s, ℓ_star; param_str="star_s", kwargs...)
		else
			star_s_σ = nothing
		end

		model.metadata[:todo][:err_estimated] = true
	    if save; @save save_fn rvs rvs_σ tel_s_σ star_s_σ end
		if save_model; @save save_model_fn model end
		return rvs, rvs_σ, tel_s_σ, star_s_σ
	else
		println("loading σs")
		if save; @load save_fn rvs rvs_σ tel_s_σ star_s_σ end
		return rvs, rvs_σ, tel_s_σ, star_s_σ
	end
end


function estimate_σ_curvature_helper(x::AbstractVecOrMat, ℓ::Function; n::Int=7, use_gradient::Bool=false, multithread::Bool=nthreads() > 3, print_every::Int=10, kwargs...)
	σs = Array{Float64}(undef, length(x))
	if !multithread
		x_test = Array{Float64}(undef, n)
		ℓs = Array{Float64}(undef, n)
	end
	if use_gradient; g = ∇(ℓ) end
	_std = std(x)
	if multithread
		nchains = nthreads()
		schedule = collect(Iterators.partition(eachindex(x), Int(ceil(length(x)/nchains))))
		# Threads.@threads for i in 1:nchains
		# ThreadsX.foreach(1:nchains) do i
		@threads for i in 1:nchains
			local _todo = copy(schedule[i])
			local _σs = Array{Float64}(undef, length(_todo))
			local _x = copy(x)
			local _x_test = Array{Float64}(undef, n)
			local _ℓs = Array{Float64}(undef, n)
			for ii in eachindex(_todo)
				k = _todo[ii]
				_x_test .= _x[k] .+ LinRange(-_std, _std, n)
				for j in 1:n
					_x[k] = _x_test[j]
					if use_gradient
						_ℓs[j] = only(g(_x))[k]
					else
						_ℓs[j] = ℓ(_x)
					end
				end
				estimate_σ_curvature_helper_finalizer!(_σs, _ℓs, _x_test, ii; use_gradient=use_gradient, print_every=100000, kwargs...)
			end
			σs[_todo] .= _σs
		end
	else
		for i in eachindex(x)
			hold = x[i]
			# x_test[:] = x[i] .+ LinRange(-_std/1e3, _std/1e3, n)
			x_test[:] = x[i] .+ LinRange(-_std, _std, n)
			for j in 1:n
				x[i] = x_test[j]
				if use_gradient
					ℓs[j] = only(g(x))[i]
				else
					ℓs[j] = ℓ(x)
				end
			end
			x[i] = hold
			# println("$i: ", ℓs .- ℓs[Int(round(n//2))])
			estimate_σ_curvature_helper_finalizer!(σs, ℓs, x_test, i; use_gradient=use_gradient, print_every=print_every, kwargs...)
		end
	end
	return reshape(σs, size(x))
end

function estimate_σ_curvature_helper_finalizer!(σs::AbstractVecOrMat, _ℓs::AbstractVector, x_test::AbstractVector, i::Int; use_gradient::Bool=false, param_str::String="", print_every::Int=10, verbose::Bool=false, show_plots::Bool=false, )
	if use_gradient
		poly_f = SSOF.ordinary_lst_sq_f(_ℓs, 1; x=x_test)
		σs[i] = sqrt(1 / poly_f.w[2])
		max_dif = maximum(abs.((poly_f.(x_test)./_ℓs) .- 1))
		if verbose; println("∇_$i: $(poly_f.w[1] + poly_f.w[2] * x[i])") end
	else
		poly_f = SSOF.ordinary_lst_sq_f(_ℓs, 2; x=x_test)
		σs[i] = sqrt(1 / (2 * poly_f.w[3]))
		max_dif = maximum(abs.((poly_f.(x_test)./_ℓs) .- 1))
		if verbose; println("∇_$i: $(poly_f.w[2] + 2 * poly_f.w[3] * x[i])") end
	end
	if show_plots
		plt = scatter(x_test, _ℓs; label="ℓ")
		plot!(x_test, poly_f.(x_test); label="polynomial fit")
		display(plt)
	end
	if max_dif > 1e-2; @warn param_str * "_σ[$i] misfit at $(round(100*max_dif; digits=2))% level" end
	if i%print_every==0; println("done with $i/$(length(σs)) " * param_str * "_σ estimates") end
	# println("done with $i/$(length(σs)) " * param_str * "_σ estimates")
end

function estimate_σ_bootstrap_reducer(shaper::AbstractArray, holder::AbstractArray, reducer::Function)
	result = Array{Float64}(undef, size(shaper, 1), size(shaper, 2))
	for i in axes(shaper, 1)
		result[i, :] .= vec(reducer(view(holder, :, i, :); dims=1))
	end
	return result
end

function estimate_σ_bootstrap_helper!(rv_holder::AbstractMatrix, tel_holder, star_holder, i::Int, mws::SSOF.ModelWorkspace, data_noise::AbstractMatrix, n::Int; verbose::Bool=true)
	time_var_tel = SSOF.is_time_variable(mws.om.tel)
	time_var_star = SSOF.is_time_variable(mws.om.star)
	_mws = typeof(mws)(copy(mws.om), copy(mws.d))
	_mws.d.flux .= mws.d.flux .+ (data_noise .* randn(size(mws.d.var)))
	improve_model!(_mws, iter=50, verbose=false)
	rv_holder[i, :] = SSOF.rvs(_mws.om)
	if time_var_tel
		tel_holder[i, :, :] .= _mws.om.tel.lm.s
	end
	if time_var_star
		star_holder[i, :, :] .= _mws.om.star.lm.s
	end
	if (verbose && i%10==0); println("done with $i/$n bootstraps") end
end

function estimate_σ_bootstrap(mws::SSOF.ModelWorkspace; recalc::Bool=false, save_fn::String="", save_model_fn::String="", n::Int=50, return_holders::Bool=false, recalc_mean::Bool=false, multithread::Bool=nthreads() > 3, verbose::Bool=true)
	save = save_fn!=""
	save_model = save_model_fn==""
	if recalc || !mws.om.metadata[:todo][:err_estimated] #|| !isfile(save_fn)
	    mws.d.var[mws.d.var.==Inf] .= 0
	    data_noise = sqrt.(mws.d.var)
	    mws.d.var[mws.d.var.==0] .= Inf

	    typeof(mws.om) <: SSOF.OrderModelWobble ?
		 	rv_holder = Array{Float64}(undef, n, length(mws.om.rv)) :
			rv_holder = Array{Float64}(undef, n, length(mws.om.rv.lm.s))

		time_var_tel = SSOF.is_time_variable(mws.om.tel)
		time_var_star = SSOF.is_time_variable(mws.om.star)
		time_var_tel ?
			tel_holder = Array{Float64}(undef, n, size(mws.om.tel.lm.s, 1), size(mws.om.tel.lm.s, 2)) :
			tel_holder = nothing
		time_var_star ?
			star_holder = Array{Float64}(undef, n, size(mws.om.star.lm.s, 1), size(mws.om.star.lm.s, 2)) :
			star_holder = nothing
		if multithread
			@threads for i in 1:n
			# # using Polyester  # same performance
			# @batch per=core for i in 1:n
			# using ThreadsX  # tiny bit better performance
			# ThreadsX.foreach(1:n) do i
				estimate_σ_bootstrap_helper!(rv_holder, tel_holder, star_holder, i, mws, data_noise, n; verbose=false)
			end
		else
			for i in 1:n
				estimate_σ_bootstrap_helper!(rv_holder, tel_holder, star_holder, i, mws, data_noise, n; verbose=verbose)
			end
	    end
		recalc_mean ? rvs = vec(mean(rv_holder; dims=1)) : rvs = SSOF.rvs(mws.om)
	    rvs_σ = vec(std(rv_holder; dims=1))

		if time_var_tel
			recalc_mean ?
				tel_s = estimate_σ_bootstrap_reducer(mws.om.tel.lm.s, tel_holder, mean) :
				tel_s = mws.om.tel.lm.s
			tel_s_σ = estimate_σ_bootstrap_reducer(mws.om.tel.lm.s, tel_holder, std)
		else
			tel_s = nothing
			tel_s_σ = nothing
		end
		if time_var_star
			recalc_mean ?
				star_s = estimate_σ_bootstrap_reducer(mws.om.star.lm.s, star_holder, mean) :
				star_s = mws.om.star.lm.s
			star_s_σ = estimate_σ_bootstrap_reducer(mws.om.star.lm.s, star_holder, std)
		else
			star_s = nothing
			star_s_σ = nothing
		end
		mws.om.metadata[:todo][:err_estimated] = true
		if save_model; @save save_model_fn model end
	    if save; @save save_fn rvs rvs_σ tel_s tel_s_σ star_s star_s_σ end
		if return_holders
			return rvs, rvs_σ, tel_s, tel_s_σ, star_s, star_s_σ, rv_holder, tel_holder, star_holder
		else
			return rvs, rvs_σ, tel_s, tel_s_σ, star_s, star_s_σ
		end
	else
		@assert isfile(save_fn)
		println("loading rvs")
		if save; @load save_fn rvs rvs_σ tel_s tel_s_σ star_s star_s_σ end
		return rvs, rvs_σ, tel_s, tel_s_σ, star_s, star_s_σ
	end
end


function how_many_comps(str::String)

	# they passed a proposed amount of parameters
	if length(str) > 0 && str[1] == '[' && str[end] == ']'
		matches = [parse(Int64, t.match) for t in eachmatch(r"-?[0-9]+", str)]
		@assert length(matches) == 2 "should only pass things of the form [n_comps_tel::Int, n_comps_star::Int]"
		n_comp_tel = matches[1]
		n_comp_star = matches[2]
		# length(matches) < 3 ? better_model = 1 : better_model = matches[3]
		# @assert 0 < better_model < 3
		use_custom_n_comp = true
		# recalc = true
		println("using $n_comp_tel telluric components and $n_comp_star stellar components")

	# using the defaults
	else
		n_comp_tel = 5
		n_comp_star = 5
		# better_model = 1  # this will be ignored later
		use_custom_n_comp = false
		println("using $n_comp_tel telluric components and $n_comp_star stellar components (the default)")
	end

	return n_comp_tel, n_comp_star, use_custom_n_comp
end

function how_many_comps(str::String, desired_order::Int)

	# they said to look for a file
	if isfile(str)
		df_n_comp = DataFrame(CSV.File(str))
		i_df = SSOF.int2ind(df_n_comp.order, desired_order)
		@assert df_n_comp[i_df, :order] == desired_order
		function look(field, default)
			try
				return df_n_comp[i_df, field]
			catch err
				if typeof(err) <: ArgumentError
					return default
				end
			end
		end
		use_custom_n_comp = look(:redo, true)
		# recalc = recalc || use_custom_n_comp
		n_comp_tel = look(:n_tel_by_eye, 5)
		n_comp_star = look(:n_star_by_eye, 5)
		# better_model = look(:better_model, 1)
		# remove_reciprocal_continuum = look(:has_reciprocal_continuum, false)
		println("using $n_comp_tel telluric components and $n_comp_star stellar components from " * str)

		return n_comp_tel, n_comp_star, use_custom_n_comp
	else
		return how_many_comps(str)
	end
end


isjld2(s::String) = (isfile(s) && s[end-4:end]==".jld2")