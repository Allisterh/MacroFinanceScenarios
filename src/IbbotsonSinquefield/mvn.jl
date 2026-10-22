# I.i.d. multivariate normal scenarios: each simulated year is an independent draw from
# N(mean, Σ), with Σ_ij = σ_i σ_j C_ij built from the std σ and correlation C in a params file.

"""
    load_mvn_params(path) -> NamedTuple

Read a multivariate normal specification from the YAML file `path`, in decimal log units:

```yaml
names: [r, rb_real, re_real]
mean: {r: 0.005, rb_real: 0.015, re_real: 0.055}
std: {r: 0.035, rb_real: 0.080, re_real: 0.173}
correlation:            # rows and columns in the order of `names`
  - [1.0, 0.64, 0.28]
  - [0.64, 1.0, 0.18]
  - [0.28, 0.18, 1.0]
start_year: 2026        # first simulated year
hist_period: [1935, 2025]  # optional: historical sample shown for comparison
```

Fields: `names`, `μ`, `σ`, `C` (correlation), `Σ` (covariance), `start_year` and
`hist_period` (`nothing` when absent). Fails if `C` is not a symmetric positive-definite
correlation matrix of matching size.
"""
function load_mvn_params(path::AbstractString)
    d = YAML.load_file(path)
    names = Symbol.(d["names"])
    n = length(names)
    μ = [Float64(d["mean"][String(v)]) for v in names]
    σ = [Float64(d["std"][String(v)]) for v in names]
    rows = d["correlation"]
    length(rows) == n && all(length(r) == n for r in rows) ||
        throw(DimensionMismatch("correlation must be $n × $n (the length of names)"))
    C = Float64[rows[i][j] for i in 1:n, j in 1:n]
    C == C' && all(==(1), diag(C)) ||
        throw(ArgumentError("correlation must be symmetric with a unit diagonal"))
    all(>(0), σ) || throw(ArgumentError("std must be positive, got $σ"))
    Σ = [σ[i] * σ[j] * C[i, j] for i in 1:n, j in 1:n]   # exactly symmetric
    isposdef(Σ) || throw(ArgumentError("correlation is not positive definite"))
    hp = get(d, "hist_period", nothing)
    hist_period = hp === nothing ? nothing : (Int(hp[1]), Int(hp[2]))
    return (; names, μ, σ, C, Σ, start_year = Int(d["start_year"]), hist_period)
end

"""
    simulate_mvn(μ, Σ, n_steps, n_scen; rng=Random.default_rng()) -> Array{Float64,3}

`n_scen` paths of `n_steps` i.i.d. draws from `N(μ, Σ)`; returns `n_vars × n_steps × n_scen`.
"""
function simulate_mvn(μ::AbstractVector{<:Real}, Σ::AbstractMatrix{<:Real}, n_steps::Integer,
                      n_scen::Integer; rng::AbstractRNG = Random.default_rng())
    n_steps >= 1 || throw(ArgumentError("n_steps must be ≥ 1, got $n_steps"))
    n_scen >= 1 || throw(ArgumentError("n_scen must be ≥ 1, got $n_scen"))
    L = cholesky(Symmetric(Matrix{Float64}(Σ))).L
    Z = randn(rng, length(μ), n_steps * n_scen)
    return reshape(μ .+ L * Z, length(μ), n_steps, n_scen)
end

"""
    simulate_mvn(p::NamedTuple, n_steps, n_scen; rng) -> Scenarios

Simulate from a specification read by `load_mvn_params`, as `Scenarios` named `p.names`
starting in `p.start_year`.
"""
simulate_mvn(p::NamedTuple, n_steps::Integer, n_scen::Integer; kwargs...) =
    Scenarios(simulate_mvn(p.μ, p.Σ, n_steps, n_scen; kwargs...), p.names, p.start_year)
