using Test
using Statistics
import Plots
using StatsBase: autocor, kurtosis, skewness
using MacroFinanceScenarios
using MacroFinanceScenarios.IbbotsonSinquefield: Scenarios
using MacroFinanceScenarios.ScenariosEvaluation
using MacroFinanceScenarios.EDA: print_table
using DataFrames: DataFrame, Not

@testset "ScenariosEvaluation" begin
    using Random: Xoshiro
    rng = Xoshiro(11)
    T, n_scen = 30, 200
    data = randn(rng, 3, T, n_scen) .* [0.02, 0.05, 0.17] .+ [0.03, 0.01, 0.07]
    data[2, :, :] .+= 0.5 .* data[1, :, :]              # correlate b with a
    sc = Scenarios(data, [:a, :b, :c], 2026)
    probs = [0.05, 0.25, 0.5, 0.75, 0.95]

    @testset "path_moments" begin
        pm = path_moments(sc, :c)
        @test pm.statistic == ["mean", "std", "skewness", "ex_kurtosis", "ar1"]
        @test propertynames(pm) == [:statistic, :mean, :p5, :p25, :p50, :p75, :p95]
        X = sc[:c]
        means = vec(mean(X; dims = 1))
        @test pm[1, :mean] ≈ mean(means)
        @test collect(pm[1, 3:end]) ≈ quantile(means, probs)
        @test pm[2, :mean] ≈ mean(std.(eachcol(X)))
        @test pm[4, :p50] ≈ median(kurtosis.(eachcol(X)))
        @test pm[5, :mean] ≈ mean(autocor(x, [1])[1] for x in eachcol(X))
        all_pm = path_moments(sc)
        @test first.(all_pm) == [:a, :b, :c]
        @test last(all_pm[3]) == pm

        hist = DataFrame(variable = ["c", "a"], mean = [0.07, 0.03], std = [0.17, 0.02],
                         skewness = [-0.5, 0.1], ex_kurtosis = [1.0, 0.0], ar1 = [0.0, 0.6],
                         min = [-0.4, 0.0])
        ph = path_moments(sc, :c; hist)
        @test propertynames(ph) == [:statistic, :historical, :mean, :p5, :p25, :p50, :p75, :p95]
        @test ph.historical == [0.07, 0.17, -0.5, 1.0, 0.0]
        @test ph[:, Not(:historical)] == pm
        @test all(ismissing, path_moments(sc, :b; hist).historical)
        @test last(path_moments(sc; hist)[1]).historical == [0.03, 0.02, 0.1, 0.0, 0.6]
    end

    @testset "path correlations" begin
        # 2 paths of 3 variables, correlations checked by hand.
        x = [1.0, 2.0, 3.0, 4.0]
        p1 = permutedims(hcat(x, 2x, -x))                    # cor(a, b) = 1, cor(a, c) = -1
        p2 = permutedims(hcat(x, [1.0, 3.0, 2.0, 4.0], x))   # cor(a, b) = 0.8, cor(a, c) = 1
        small = Scenarios(cat(p1, p2; dims = 3), [:a, :b, :c], 2026)
        t = path_correlations(small; probs = [0.0, 1.0])
        @test propertynames(t) == [:var1, :var2, :mean, :p0, :p100]
        @test t.var1 == ["a", "a", "b"] && t.var2 == ["b", "c", "c"]   # lower triangle, by column
        @test t.mean ≈ [0.9, 0.0, -0.1]                                # cor(b, c): -1 and 0.8
        @test t.p0 ≈ [0.8, -1.0, -1.0] && t.p100 ≈ [1.0, 1.0, 0.8]

        hist = DataFrame(variable = ["c", "a"], c = [1.0, 0.25], a = [0.25, 1.0])
        th = path_correlations(small; hist, probs = [0.0, 1.0])
        @test propertynames(th) == [:var1, :var2, :historical, :mean, :p0, :p100, :hist_rank]
        @test isequal(th.historical, [missing, 0.25, missing])        # b not in hist
        @test isequal(th.hist_rank, [missing, 0.5, missing])          # paths at -1 and 1
        @test th[:, Not([:historical, :hist_rank])] == t

        m = mean_path_correlations(small)
        @test propertynames(m) == [:variable, :a, :b, :c]
        @test m.variable == ["a", "b", "c"]
        M = Matrix(m[:, 2:end])
        @test M ≈ M' && all(M[i, i] ≈ 1 for i in 1:3)
        @test M[2, 1] ≈ 0.9 && M[3, 1] ≈ 0.0 && M[3, 2] ≈ -0.1

        # Random scenarios: default percentiles and the matrix agree with the table's mean.
        pc = path_correlations(sc)
        @test propertynames(pc) == [:var1, :var2, :mean, :p5, :p25, :p50, :p75, :p95]
        @test pc.mean[1] ≈ mean(cor(permutedims(data[:, :, s]))[2, 1] for s in 1:n_scen)
        @test mean_path_correlations(sc).a[2:3] ≈ pc.mean[1:2]
    end

    @testset "horizon values and percentiles" begin
        X = sc[:a]
        @test horizon_values(X, 5) ≈ vec(sum(X[1:5, :]; dims = 1)) ./ 5
        @test horizon_values(X, 5; kind = :level) == X[5, :]
        @test horizon_values(X, 5; kind = :wealth) ≈ exp.(vec(sum(X[1:5, :]; dims = 1)))
        @test horizon_values(X, 1) ≈ X[1, :]
        @test_throws ArgumentError horizon_values(X, T + 1)
        @test_throws ArgumentError horizon_values(X, 5; kind = :foo)
        # simple = true: exp(·) − 1 per path, before any summary.
        @test horizon_values(X, 5; simple = true) ≈ expm1.(horizon_values(X, 5))
        @test horizon_values(X, 5; kind = :level, simple = true) ≈ expm1.(X[5, :])
        @test_throws ArgumentError horizon_values(X, 5; kind = :wealth, simple = true)
        hs = horizon_percentiles(sc; vars = [:c], horizons = [10], simple = true)
        hl = horizon_percentiles(sc; vars = [:c], horizons = [10])
        @test hs.mean[1] ≈ mean(expm1.(horizon_values(sc[:c], 10)))
        @test hs.mean[1] > expm1(hl.mean[1])                 # Jensen: mean does not commute
        # Percentiles commute (up to the interpolation between order statistics).
        @test hs.p50[1] ≈ expm1(hl.p50[1]) rtol = 1e-4
        ext = horizon_percentiles(sc; vars = [:c], horizons = [10], probs = [0.0, 1.0], simple = true)
        @test ext.p0[1] ≈ expm1(minimum(horizon_values(sc[:c], 10)))

        hp = horizon_percentiles(sc; vars = [:a, :c], horizons = [1, 10, 25])
        @test hp.variable == ["a", "a", "a", "c", "c", "c"]
        @test hp.horizon == [1, 10, 25, 1, 10, 25]
        @test hp.year == [2026, 2035, 2050, 2026, 2035, 2050]
        @test collect(hp[5, [:p5, :p25, :p50, :p75, :p95]]) ≈ quantile(horizon_values(sc[:c], 10), probs)
        @test hp[5, :mean] ≈ mean(horizon_values(sc[:c], 10))
        # Annualised dispersion shrinks with the horizon for i.i.d. returns.
        @test hp[6, :p95] - hp[6, :p5] < hp[4, :p95] - hp[4, :p5]

        fan = horizon_percentiles(sc; vars = [:a], horizons = 1:T, kind = :level, probs = [0.025, 0.5])
        @test propertynames(fan) == [:variable, :horizon, :year, :mean, Symbol("p2.5"), :p50]
        @test fan.p50 ≈ [median(sc[:a][t, :]) for t in 1:T]
    end

    @testset "fan" begin
        @test fan(sc, :a) isa Plots.Plot
        @test fan(sc, :a; kind = :annualised, simple = true, μ = 0.03, scale = 100) isa Plots.Plot
        p = fan(sc, :c; kind = :wealth, μ = 0.07, yscale = :log10, title = "c")
        @test p isa Plots.Plot
        @test last(p.series_list)[:y] ≈ exp.((1:T) .* 0.07)  # reference curve exp(h μ)
    end

    @testset "drawdowns" begin
        r = log.([1.1 0.9 1.0; 0.8 1.0 1.0; 1.2 1.0 1.0; 1.0 1.2 1.0])  # wealth paths per column
        dd, len = max_drawdown_and_length(r)
        # col 1: 1.1 → 0.88 → 1.056 → 1.056, peak 1.1, under water 3 years
        @test dd[1] ≈ 0.2
        @test len[1] == 3
        # col 2: loss in year 1 counts (wealth starts at 1); 0.9 ×1 ×1 ×1.2 = 1.08 recovers
        @test dd[2] ≈ 0.1
        @test len[2] == 3
        @test dd[3] == 0 && len[3] == 0

        dt = drawdown_table(sc; vars = [:c])
        @test dt.statistic == ["max_drawdown", "max_dd_length"]
        d, l = max_drawdown_and_length(sc[:c])
        @test dt[1, :mean] ≈ mean(d)
        @test dt[2, :p50] ≈ median(l)
        @test all(0 .<= d .< 1)
    end

    @testset "annualise and period_percentiles" begin
        X = reshape(1.0:12.0, 2, 6)                     # rows = scenarios, cols = periods
        @test annualise(X, 2) == [4.0 12.0 20.0; 6.0 14.0 22.0]
        @test annualise(X, 4) == reshape([16.0, 20.0], 2, 1)
        pp = period_percentiles(X, [0.0, 1.0]; freq = 2)
        @test pp.period == [1, 2, 3]
        @test pp.p0 == [4.0, 12.0, 20.0] && pp.p100 == [6.0, 14.0, 22.0]
    end

    @test occursin("ex_kurtosis", sprint(io -> print_table(io, path_moments(sc, :a))))
end
