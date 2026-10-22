# AR(1) fitted by OLS: x_t = c + φ x_{t-1} + e_t. Dependency-free (plain least squares).

"""
    AR1

OLS fit of `x_t = c + φ x_{t-1} + e_t`.

- `c`, `φ`    coefficients; `se = (se_c, se_φ)` their OLS standard errors
- `σ`         residual std, with `n − 2` degrees of freedom
- `r2`        R² of the regression
- `resid`, `fitted`, `years`: one entry per observation `t = 2..n` (the first is lost to the lag)
- `period`    `(first_year, last_year)` of the sample passed to `fit_ar1`, including the lag year
"""
struct AR1
    c::Float64
    φ::Float64
    σ::Float64
    se::NTuple{2,Float64}
    r2::Float64
    resid::Vector{Float64}
    fitted::Vector{Float64}
    years::Vector{Int}
    period::Tuple{Int,Int}
end

"""
    fit_ar1(x::AbstractVector; years=1:length(x)) -> AR1
    fit_ar1(ta::TimeArray) -> AR1
    fit_ar1(ta::TimeArray, var::Symbol) -> AR1

Fit an AR(1) by OLS. With a `TimeArray` (one column, or column `var`), the years come from
its timestamps, so calling it on a `subperiod` slice records that period.
"""
function fit_ar1(x::AbstractVector{<:Real}; years::AbstractVector{<:Integer} = 1:length(x))
    n = length(x)
    n >= 3 || throw(ArgumentError("fit_ar1 needs at least 3 observations, got $n"))
    length(years) == n || throw(DimensionMismatch("years has length $(length(years)), x has length $n"))

    y = x[2:end]
    X = [ones(n - 1) x[1:end-1]]
    β = X \ y
    fitted = X * β
    resid = y .- fitted

    σ² = sum(abs2, resid) / (n - 1 - 2)
    se = sqrt.(σ² .* diag(inv(X' * X)))
    r2 = 1 - sum(abs2, resid) / sum(abs2, y .- mean(y))

    return AR1(β[1], β[2], sqrt(σ²), (se[1], se[2]), r2, resid, fitted,
               collect(Int, years[2:end]), (Int(first(years)), Int(last(years))))
end

function fit_ar1(ta::TimeArray)
    length(colnames(ta)) == 1 || throw(ArgumentError("fit_ar1 needs a single-column TimeArray; pass the column name"))
    return fit_ar1(vec(values(ta)); years = year.(timestamp(ta)))
end

fit_ar1(ta::TimeArray, var::Symbol) = fit_ar1(ta[var])

"""Unconditional mean `c / (1 − φ)` of a stationary AR(1)."""
unconditional_mean(ar::AR1) = ar.c / (1 - ar.φ)

"""Unconditional std `σ / √(1 − φ²)` of a stationary AR(1)."""
unconditional_std(ar::AR1) = ar.σ / sqrt(1 - ar.φ^2)

"""
    simulate_ar1(ar::AR1, T, n_scen; y0=:mean, rng=Random.default_rng()) -> Matrix{Float64}

`T × n_scen` paths of `x_t = ar.c + ar.φ x_{t-1} + ar.σ ε_t` with i.i.d. `ε_t ~ N(0, 1)`.
`y0` is the value before the first simulated period: `:mean` (the unconditional mean) or a number.
"""
function simulate_ar1(ar::AR1, T::Integer, n_scen::Integer;
                      y0::Union{Symbol,Real} = :mean, rng::AbstractRNG = Random.default_rng())
    sc = simulate_ar1([ar], [:x], T, n_scen; y0 = y0 isa Symbol ? y0 : [y0], rng)
    return sc.data[1, :, :]
end

"""
    simulate_ar1(models::AbstractVector{AR1}, names, T, n_scen;
                 C=I, y0=:mean, start_year=1, rng=Random.default_rng()) -> Scenarios

Simulate `n_scen` paths of length `T` from independent AR(1) equations
`x_{i,t} = c_i + φ_i x_{i,t-1} + σ_i ε_{i,t}`, one per model, rows labelled by `names`. The
shocks `ε_t` are i.i.d. `N(0, C)` across years, with `C` the correlation matrix of the
standardised residuals (identity by default). `y0` is the state before the first simulated
period: `:mean` (the unconditional means) or a vector.

Throws if some `|φ_i| ≥ 1`.
"""
function simulate_ar1(models::AbstractVector{AR1}, names::AbstractVector{Symbol}, T::Integer,
                      n_scen::Integer; C::AbstractMatrix{<:Real} = Matrix{Float64}(I, length(models), length(models)),
                      y0::Union{Symbol,AbstractVector{<:Real}} = :mean, start_year::Integer = 1,
                      rng::AbstractRNG = Random.default_rng())
    n = length(models)
    length(names) == n || throw(DimensionMismatch("$n models, names has $(length(names))"))
    size(C) == (n, n) || throw(DimensionMismatch("C must be $n × $n, got $(size(C))"))
    T >= 1 || throw(ArgumentError("T must be ≥ 1, got $T"))
    n_scen >= 1 || throw(ArgumentError("n_scen must be ≥ 1, got $n_scen"))
    all(m -> abs(m.φ) < 1, models) ||
        throw(ArgumentError("every AR(1) must be stationary, got φ = $([m.φ for m in models])"))

    c = [m.c for m in models]
    φ = [m.φ for m in models]
    σ = [m.σ for m in models]
    y_init = if y0 === :mean
        unconditional_mean.(models)
    elseif y0 isa Symbol
        throw(ArgumentError("y0 must be :mean or a vector, got :$y0"))
    else
        length(y0) == n || throw(DimensionMismatch("y0 has length $(length(y0)), expected $n"))
        Vector{Float64}(y0)
    end

    L = cholesky(Symmetric(Matrix{Float64}(C))).L
    sims = reshape(σ .* (L * randn(rng, n, T * n_scen)), n, T, n_scen)   # shocks, then levels in place
    for s in 1:n_scen
        @views sims[:, 1, s] .+= c .+ φ .* y_init
        for t in 2:T
            @views sims[:, t, s] .+= c .+ φ .* sims[:, t-1, s]
        end
    end
    return Scenarios(sims, names, start_year)
end
