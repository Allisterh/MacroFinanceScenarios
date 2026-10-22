# Fit the AR(1) scenario model on the annual Shiller data and save it as YAML (`save_ar1_params`)
# for ar1_report.qmd. State: the real short rate `r` (AR(1) fit), the excess bond return `tp`
# and the excess equity return `xr` (intercept only: φ = 0, c = the sample mean). Shock
# correlations come from the residuals over the same years; `y0` is the last observation.
#
# Usage (from the repo root):
#   julia --project analisys/Ibbotson-Sinquefield-SimulationModel/export_ar1_params.jl [y0 y1]

using CSV
using Dates: year
using Statistics
using TimeSeries
using MacroFinanceScenarios.IbbotsonSinquefield

const ROOT = joinpath(@__DIR__, "..", "..")

period = length(ARGS) >= 2 ? (parse(Int, ARGS[1]), parse(Int, ARGS[2])) : (1935, 2025)
data = subperiod(TimeArray(CSV.File(joinpath(ROOT, "data", "shiller_annual.csv")); timestamp = :timestamp),
                 period)

# Intercept-only AR(1) on the same years as the lagged fit (the first year is dropped), so
# the residuals line up and every equation has the same n.
function fit_mean(ta::TimeArray)
    yrs = year.(timestamp(ta))
    x = vec(values(ta))[2:end]
    n, c, σ = length(x), mean(x), std(x)
    return AR1(c, 0.0, σ, (σ / sqrt(n), 0.0), 0.0, x .- c, fill(c, n), yrs[2:end], (yrs[1], yrs[end]))
end

names = [:r, :tp, :xr]
models = [fit_ar1(data, :r), fit_mean(data[:tp]), fit_mean(data[:xr])]
C = cor(stack(m.resid for m in models))
y0 = vec(values(data[names...])[end, :])

out = joinpath(@__DIR__, "specs", "ar1", "hist$(period[1])-$(period[2]).yaml")
save_ar1_params(out; models, names, C, y0, start_year = period[2] + 1,
                header = """
                AR(1) scenario model, decimal log units (export_ar1_params.jl, save_ar1_params).
                r: real T-bill return, AR(1) fit; tp: excess bond return rb − rf; xr: excess equity
                return re − rf; tp and xr are intercept only (phi = 0, c = sample mean).""")
println("Wrote $out")
