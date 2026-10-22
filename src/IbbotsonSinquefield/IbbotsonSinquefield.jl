module IbbotsonSinquefield

export load_shiller_annual, subperiod
export ColumnTransformer, fit_columns, transform_columns, inverse_transform_columns
export Except, Affine, Invertible, Demean, Standardize
export AR1, fit_ar1, unconditional_mean, unconditional_std, simulate_ar1
export Scenarios, scenario_years, is_coefficients, issm_bootstrap_matrix
export simulate_var_bootstrapped_res, returns_from_components
export save_issm_params, load_issm_params, save_residuals, load_residuals, save_ar1_params, load_ar1_params
export match_means, block_bootstrap
export load_mvn_params, simulate_mvn

using CSV
using Dates
using Distributions: Geometric
using LinearAlgebra: Diagonal, I, Symmetric, cholesky, diag, eigvals, isposdef, mul!
using Random
using Statistics
using TimeSeries
using XLSX
using YAML

include(joinpath(@__DIR__, "..", "common", "bond_returns.jl"))
include(joinpath(@__DIR__, "data.jl"))
include(joinpath(@__DIR__, "ar1.jl"))
include(joinpath(@__DIR__, "simulate.jl"))
include(joinpath(@__DIR__, "transform.jl"))
include(joinpath(@__DIR__, "moments.jl"))
include(joinpath(@__DIR__, "block_bootstrap.jl"))
include(joinpath(@__DIR__, "mvn.jl"))
include(joinpath(@__DIR__, "io.jl"))

end
