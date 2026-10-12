# Scenario evaluation: per-path moments and correlations, horizon percentiles, fan charts and
# drawdowns of simulated `Scenarios` (n_vars × T × n_scen, log units). Like the EDA
# helpers, every table function returns a labelled DataFrame and never prints; use
# `EDA.print_table` for the output.

const DEFAULT_PROBS = [0.05, 0.25, 0.5, 0.75, 0.95]

# Column name for probability p: 0.05 → :p5, 0.025 → :p2.5.
function prob_name(p::Real)
    x = round(100p; digits = 6)
    return Symbol("p", isinteger(x) ? Int(x) : x)
end

# NamedTuple-ready pairs: mean of `x` followed by its quantiles at `probs`.
summary_pairs(x::AbstractVector, probs) =
    [:mean => mean(x); [prob_name(p) => q for (p, q) in zip(probs, quantile(x, probs))]]

"""
    path_moments(sc::Scenarios, name; hist=nothing, probs=[0.05, 0.25, 0.5, 0.75, 0.95]) -> DataFrame
    path_moments(sc::Scenarios; hist=nothing, probs) -> Vector{Pair{Symbol,DataFrame}}

Moments of variable `name` computed within each simulated path (mean, std, skewness,
excess kurtosis, lag-1 autocorrelation, `neg`: percent of negative values), summarised across scenarios: one row per
statistic, columns `mean` and the percentiles `probs` of that statistic. Without `name`,
one table per variable of `sc`, in order.

`hist` is an optional `EDA.describe_series` table (in the same units as `sc`): its row
with `variable == name` goes into a column `historical` after `statistic`, `missing` if
`hist` has no such row.
"""
function path_moments(sc::Scenarios, name::Symbol; hist = nothing, probs = DEFAULT_PROBS)
    X = sc[name]                                        # T × n_scen
    size(X, 1) >= 2 || throw(ArgumentError("path moments need at least 2 simulated years"))
    paths = eachcol(X)
    stats = [:mean => mean.(paths),
             :std => std.(paths),
             :skewness => skewness.(paths),
             :ex_kurtosis => kurtosis.(paths),          # StatsBase.kurtosis is excess kurtosis
             :ar1 => [autocor(x, [1])[1] for x in paths],
             :neg => [100 * mean(<(0), x) for x in paths]]
    rows = [(; :statistic => String(s), summary_pairs(v, probs)...) for (s, v) in stats]
    df = DataFrame(rows)
    if hist !== nothing
        i = findfirst(==(String(name)), String.(hist.variable))
        insertcols!(df, 2, :historical =>
            [i === nothing ? missing : hist[i, s] for s in first.(stats)])
    end
    return df
end

path_moments(sc::Scenarios; hist = nothing, probs = DEFAULT_PROBS) =
    [name => path_moments(sc, name; hist, probs) for name in sc.names]

# n × n × n_scen correlations of the variables of `sc`, computed within each path over time.
function path_correlation_array(sc::Scenarios)
    n, T, n_scen = size(sc.data)
    T >= 2 || throw(ArgumentError("correlations need at least 2 simulated years"))
    C = Array{Float64,3}(undef, n, n, n_scen)
    for s in 1:n_scen
        C[:, :, s] = cor(permutedims(@view sc.data[:, :, s]))
    end
    return C
end

"""
    path_correlations(sc::Scenarios; hist=nothing, probs=[0.05, 0.25, 0.5, 0.75, 0.95]) -> DataFrame

Correlations of each pair of variables of `sc` computed within each path (over time),
summarised across scenarios: one row per pair `var1`, `var2` in lower-triangle order
(column by column of the correlation matrix), columns `mean` and the percentiles `probs`.

`hist` is an optional `EDA.correlation_table` (matched by variable name): it adds a column
`historical` after `var2` and a last column `hist_rank`, the share of paths with a
correlation at or below the historical one; both `missing` for a pair not in `hist`.
"""
function path_correlations(sc::Scenarios; hist = nothing, probs = DEFAULT_PROBS)
    C = path_correlation_array(sc)
    n = length(sc.names)
    labels = String.(sc.names)
    pairs = [(i, j) for j in 1:n for i in j+1:n]
    df = DataFrame([(; :var1 => labels[j], :var2 => labels[i], summary_pairs(C[i, j, :], probs)...)
                    for (i, j) in pairs])
    if hist !== nothing
        function hist_value(a, b)
            k = findfirst(==(a), String.(hist.variable))
            return k === nothing || !(b in names(hist)) ? missing : hist[k, b]
        end
        h = [hist_value(labels[j], labels[i]) for (i, j) in pairs]
        insertcols!(df, 3, :historical => h)
        df.hist_rank = [ismissing(x) ? missing : mean(<=(x), C[i, j, :])
                        for (x, (i, j)) in zip(h, pairs)]
    end
    return df
end

"""
    mean_path_correlations(sc::Scenarios) -> DataFrame

Mean across paths of the within-path correlation matrices of `sc` (the `mean` column of
`path_correlations` as a matrix), in the layout of `EDA.correlation_table`: column
`variable` with the row labels, then one column per variable.
"""
function mean_path_correlations(sc::Scenarios)
    labels = String.(sc.names)
    df = DataFrame(dropdims(mean(path_correlation_array(sc); dims = 3); dims = 3), labels)
    insertcols!(df, 1, :variable => labels)
    return df
end

"""
    horizon_values(X::AbstractMatrix, h; kind=:annualised, simple=false) -> Vector

Value at horizon `h` of each path (column) of the `T × n_scen` log series `X`:

- `:annualised`: annualised cumulative `sum(x_1..h) / h` (returns, average inflation)
- `:level`: the level in year h, `x_h` (state variables: π, r, rf, tp)
- `:wealth`: cumulative wealth `exp(sum(x_1..h))` (or the price level for π)

`simple = true` turns `:level` and `:annualised` log rates into simple rates, `exp(·) − 1`
(the geometric-average annual return for `:annualised`). It is applied per path, so summary
statistics of the result are those of the simple rates (the mean does not commute with
`exp`). Not allowed with `:wealth`, which is already exponentiated.
"""
function horizon_values(X::AbstractMatrix, h::Integer; kind::Symbol = :annualised, simple::Bool = false)
    1 <= h <= size(X, 1) || throw(ArgumentError("horizon $h outside 1:$(size(X, 1))"))
    kind in (:annualised, :level, :wealth) ||
        throw(ArgumentError("kind must be :annualised, :level or :wealth, got :$kind"))
    kind === :wealth && simple && throw(ArgumentError("simple = true does not apply to kind = :wealth"))
    v = if kind === :level
        X[h, :]
    else
        S = vec(sum(@view(X[1:h, :]); dims = 1))
        kind === :wealth ? exp.(S) : S ./ h
    end
    return simple ? expm1.(v) : v
end

"""
    horizon_percentiles(sc::Scenarios; vars=sc.names, horizons=[1, 5, 10, 25],
                        probs=[0.05, 0.25, 0.5, 0.75, 0.95], kind=:annualised,
                        simple=false) -> DataFrame

Distribution across scenarios of each variable in `vars` at each horizon (see
`horizon_values` for `kind` and `simple`): one row per variable × horizon with columns
`variable`, `horizon`, `year`, `mean` and the percentiles `probs`. Horizons beyond the
simulated length throw.

Which `kind` fits which variable: `:annualised` for asset returns (rf, rb, re), π, erp;
`:level` for the persistent state variables π, r, rf and tp (with `horizons = 1:T` this is
the fan chart data).
"""
function horizon_percentiles(sc::Scenarios; vars = sc.names, horizons = [1, 5, 10, 25],
                             probs = DEFAULT_PROBS, kind::Symbol = :annualised, simple::Bool = false)
    rows = map(Iterators.product(horizons, vars)) do (h, name)
        v = horizon_values(sc[name], h; kind, simple)
        (; :variable => String(name), :horizon => h, :year => sc.start_year + h - 1,
           summary_pairs(v, probs)...)
    end
    return DataFrame(vec(rows))
end

"""
    fan(sc::Scenarios, name; kind=:level, simple=false, μ=nothing, scale=1.0,
        title=String(name), yscale=:identity) -> Plots.Plot

Fan chart of variable `name` over the simulated years: 5–95 and 25–75 percentile bands and
the median of `horizon_values(sc[name], h; kind, simple)` for every `h`, multiplied by
`scale` (e.g. 100 for percent; keep 1 for `:wealth`).

`μ` is an optional reference line in log units (e.g. the unconditional mean), drawn dashed
at `μ` for `:level`/`:annualised`, `exp(μ) − 1` with `simple = true`, and along `exp(h μ)`
for `:wealth`. `yscale` is passed to Plots, e.g. `:log10` for long wealth fans.
"""
function fan(sc::Scenarios, name::Symbol; kind::Symbol = :level, simple::Bool = false,
             μ::Union{Real,Nothing} = nothing, scale::Real = 1.0, title = String(name),
             yscale::Symbol = :identity)
    q = horizon_percentiles(sc; vars = [name], horizons = 1:size(sc.data, 2), kind, simple)
    yrs = q.year
    p = plot(yrs, scale .* q.p5; fillrange = scale .* q.p95, fillalpha = 0.2, lw = 0, color = 1,
             label = "5–95", title, legend = :topleft, legendfontsize = 7, yscale)
    plot!(p, yrs, scale .* q.p25; fillrange = scale .* q.p75, fillalpha = 0.3, lw = 0, color = 1,
          label = "25–75")
    plot!(p, yrs, scale .* q.p50; lw = 2, color = 1, label = "median")
    if μ !== nothing
        ref = kind === :wealth ? exp.(q.horizon .* μ) : fill(simple ? expm1(μ) : float(μ), length(yrs))
        plot!(p, yrs, scale .* ref; ls = :dash, color = :black, label = "uncond. mean")
    end
    return p
end

"""
    max_drawdown_and_length(returns::AbstractMatrix) -> (max_drawdowns, max_dd_lengths)

Per path (column) of the `T × n_scen` log returns: the maximum drawdown of wealth
`exp(cumsum(returns))` from its running peak (a fraction, 0 = never below the peak) and the
longest drawdown, the number of consecutive years spent below the previous peak. Wealth
starts at 1 before the first year, so a loss in year 1 already counts as a drawdown.
"""
function max_drawdown_and_length(returns::AbstractMatrix)
    T, n_scen = size(returns)
    max_dd = zeros(n_scen)
    max_len = zeros(Int, n_scen)
    for s in 1:n_scen
        logw, peak, len = 0.0, 0.0, 0          # log wealth and its running peak
        for t in 1:T
            logw += returns[t, s]
            if logw >= peak
                peak, len = logw, 0
            else
                len += 1
                max_dd[s] = max(max_dd[s], 1 - exp(logw - peak))
                max_len[s] = max(max_len[s], len)
            end
        end
    end
    return max_dd, max_len
end

"""
    drawdown_table(sc::Scenarios; vars=sc.names, probs=[0.05, 0.25, 0.5, 0.75, 0.95]) -> DataFrame

Maximum drawdown and longest drawdown length (`max_drawdown_and_length`) per path for each
return variable in `vars`, summarised across scenarios: rows `variable` × `statistic`
(`max_drawdown`, `max_dd_length`), columns `mean` and the percentiles `probs`. Pass nominal
and real returns (`returns_from_components`) separately.
"""
function drawdown_table(sc::Scenarios; vars = sc.names, probs = DEFAULT_PROBS)
    rows = NamedTuple[]
    for name in vars
        dd, len = max_drawdown_and_length(sc[name])
        for (stat, v) in (:max_drawdown => dd, :max_dd_length => float.(len))
            push!(rows, (; :variable => String(name), :statistic => String(stat),
                         summary_pairs(v, probs)...))
        end
    end
    return DataFrame(rows)
end

"""
    annualise(scenarios::AbstractMatrix, shift=2) -> Matrix

Sum non-overlapping blocks of `shift` columns of `scenarios` (rows = scenarios, columns =
periods), e.g. quarterly log returns into annual ones with `shift = 4`. A trailing
incomplete block is dropped.
"""
function annualise(scenarios::AbstractMatrix, shift::Integer = 2)
    n_blocks = size(scenarios, 2) ÷ shift
    result = zeros(size(scenarios, 1), n_blocks)
    for p in 1:n_blocks
        result[:, p] .= vec(sum(@view(scenarios[:, (p-1)*shift+1:p*shift]); dims = 2))
    end
    return result
end

"""
    period_percentiles(X::AbstractMatrix, probs; freq=1) -> DataFrame

Percentiles `probs` across scenarios of each period of `X` (rows = scenarios, columns =
periods), after summing blocks of `freq` periods (`annualise`). One row per aggregated
period. Replaces TCVAR's `print_percentiles`: it returns the table instead of printing it.
"""
function period_percentiles(X::AbstractMatrix, probs; freq::Integer = 1)
    Y = annualise(X, freq)
    rows = [(; :period => p, [prob_name(q) => v for (q, v) in zip(probs, quantile(y, probs))]...)
            for (p, y) in enumerate(eachcol(Y))]
    return DataFrame(rows)
end
