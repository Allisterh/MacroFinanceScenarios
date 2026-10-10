### Refactor ISSM workflow and evaluation report

Refactor [issm_report.qmd](issm_report.qmd), [issm_eda.ipynb](issm_eda.ipynb),
[render_reports.sh](render_reports.sh), [IbbotsonSinquefield](../../src/IbbotsonSinquefield) and
[ScenariosEvaluation](../../src/ScenariosEvaluation/evaluation.jl) to the workflow below.

## Target workflow

```
issm_eda.ipynb  ──estimates, writes──▶  specs/<spec>/ ──read by──▶  issm_report.qmd  ──▶  reports/<spec>.md
                                          params.yaml               (no fitting, no ie_data.xlsx)
                                          bootstrap_residuals.csv
                                          ar_residuals.csv
                                          hist_moments.csv
                                          hist_correlations.csv
```

- **issm_eda** is the only place where data is loaded and models are estimated. A final
  *Export* section writes one spec directory per specification
  (`specs/ar1935-2025_boot1935-2025_erp4.0/`, same naming as the current reports).
- **issm_report** only loads a spec directory, simulates and presents. Its Quarto params
  shrink to `spec` (directory name) plus the simulation settings `T`, `n_scen`, `seed`.
- **render_reports.sh** loops over `specs/*/` instead of the hard-coded `specs=(...)` list.
- **Units:** everything saved to disk and everything passed to library functions is in
  decimal log units. ×100 happens only when building a table or plot.

## 1. Model parameters in YAML (estimated in EDA, presented in report)

`specs/<spec>/params.yaml`, all values in decimals:

```yaml
ar_period: [1935, 2025]
boot_period: [1935, 2025]
erp_target: 0.04
last_year: 2025                 # simulation starts in last_year + 1
variables: [π, r, tp, erp]
y0: [0.027, 0.015, 0.004, 0.08] # initial state = last observation of [π, r, tp, erp]
c:   [...]                      # VAR intercepts (is_coefficients)
phi: [..., ..., ..., 0.0]       # diagonal of A
ar1:                            # for the presentation table only
  π:  {n: 90, c: ..., se_c: ..., phi: ..., se_phi: ..., sigma: ..., r2: ...}
  r:  {...}
  tp: {...}
```

- Add `save_issm_params(path; ...)` / `load_issm_params(path) -> (; c, A, names, y0, start_year, ar_table, ...)`
  to `IbbotsonSinquefield` (new `io.jl`), so the schema lives in one place and has a
  round-trip test. New dependency: `YAML.jl`.
- Report: replace the `fit_ar1` / `is_coefficients` calls with `load_issm_params`; call
  `simulate_var_bootstrapped_res(...; y0 = p.y0, start_year = p.start_year)` (no `data`).
  The *AR(1) fits* and *VAR(1) intercepts* tables are built from the loaded values.

## 2. Residuals (estimated in EDA, saved descaled)

EDA scales data ×100 at load. For the export, **refit on the unscaled `raw` data** rather
than dividing percent results by 100 (exact, and no risk of descaling φ or R² by mistake).

- `bootstrap_residuals.csv`: `year` + columns `π, r, tp, erp` = the demeaned VAR residuals
  `U` from `issm_bootstrap_matrix` over the bootstrap window. This is what the simulation
  draws from.
- `ar_residuals.csv`: `year` + `π, r, tp` = raw AR(1) residuals (`fit.resid`, not
  demeaned) over the AR window. For diagnostics/presentation only.

Two files because the AR and bootstrap windows can differ (different year ranges).
Format: CSV via `CSV.jl` (human-readable, round-trips through DataFrames; new dependency).

## 3. Real bond and equity returns in `data.jl`

- Add `rb_real = rb − π` and `re_real = re − π` to `load_shiller_annual` (docstring +
  test `v(:rb) .- v(:π) ≈ v(:rb_real)` etc.). The real T-bill return is already `r`.
- Rename the outputs of `returns_from_components(sc; real = true)` from `[:rf, :rb, :re]`
  to `[:r, :rb_real, :re_real]`, so simulated and historical series share names (needed
  for 6 and 7). Update its test.
- EDA: add `rb_real`, `re_real` to `eda_vars` (moments and correlations per period).

## 4. Historical moments (saved in EDA)

Over the **bootstrap window** of the spec (the sample the simulation reproduces), in
decimals:

- `hist_moments.csv`: `describe_series` output for `π, r, rf, tp, rb, re, rb_real,
  re_real, erp` (historical series) and for the columns of `U` with names prefixed `u_`
  (`u_π, u_r, u_tp, u_erp`).
- `hist_correlations.csv`: `correlation_table` of the same historical series.

The report's *Historical moments* section and the *Validation* section use these files
instead of `describe_series(boot_data)` / `cor(values(boot_data))`.

## 5. Scaling only at presentation

Keep scenarios and historical tables in decimals; scale with the existing `pct` helper
(for `Scenarios`) plus a new `pct(df::DataFrame; cols)` for the loaded moment tables.
Means, std and percentiles ×100; skewness, kurtosis, autocorrelations, correlations not.
(The report already does this for scenarios; the change is to apply it to the loaded
historical tables too and never to data used in computation.)

**Annualised values and levels are presented as simple rates**, not log rates:

- annualised: `exp(sum(x_1..h)/h) − 1` (geometric-average annual return, or annual
  inflation rate for `π`)
- level in year h: `exp(x_h) − 1` (simple inflation, short rates, term premium)

Add a `simple::Bool = false` kwarg to `horizon_values` / `horizon_percentiles` that
applies `exp(·) − 1` for `kind = :level` and `:annualised` (throws for `:wealth`, which is
already exponentiated). It converts **per path before summarising**: percentiles commute
with the monotone `exp(·) − 1`, but the mean does not (`mean(exp(x)) ≠ exp(mean(x))`), so
converting a finished log table is wrong. The computation stays in logs; `simple = false`
is the default for library use and tests.

## 6. `fan` → `ScenariosEvaluation`

```julia
fan(sc::Scenarios, name; kind = :level, simple = false, μ = nothing, scale = 1.0,
    title = String(name)) -> Plots.Plot
```

- `kind` and `simple` as in `horizon_values` (`:level`, `:annualised`, `:wealth`), so
  the same function draws all fans in 10. Cumulative fans use the existing `:wealth`
  (`exp(sum(x_1..h))`, wealth of 1 invested).
- `μ` optional reference line, given in log units: drawn at `μ` for
  `:level`/`:annualised`, `exp(μ) − 1` with `simple = true`, and `exp(h·μ)` for
  `:wealth`. `scale` replaces the report-local `pct`; wealth fans use
  `scale = 1` (wealth is not a percentage).
- Wealth axis linear by default; optional `yscale = :log` kwarg passed through to Plots
  for long horizons.
- Adds `Plots` to the `ScenariosEvaluation` module imports (already a package dependency).
  Smoke test: returns a plot without throwing.

## 7. `path_moments` with historical column

```julia
path_moments(sc, name; hist = nothing, probs)   # hist: describe_series DataFrame
```

If `hist` is given, insert column `historical` after `statistic`, taken from the row of
`hist` with `variable == name` (statistics `mean, std, skewness, ex_kurtosis, ar1` match
`describe_series` column names). Missing variable → `missing`, not an error (e.g. `erp`
uses `erp_target` by design, still present).

## 8. `mean_correlations` → `scenarios_correlations`

```julia
scenarios_correlations(sc; hist = nothing, probs = (0.05, 0.95), digits = 2) -> DataFrame
```

- Per-path correlation matrices (`n × n × n_scen`).
- Lower triangle: `"0.30 (0.12, 0.44)"` = mean across paths and the `probs` percentiles.
- Upper triangle: historical correlation from `hist` (a `correlation_table` DataFrame,
  matched by variable name), `""` if `hist` is `nothing`.
- Diagonal: `"1"`.
- Cells are strings, so formatting happens here (`digits`), not in `print_table`.
  Update the export list, docstring and tests (ordering, bracket values on a small
  hand-checked example).

## 9. Horizon percentiles: one table per variable

Rows = horizons `[1, 5, 10, 20, 30] ∩ 1:T`, columns `year, mean, p5 … p95`, ×100.

| table | variables | `kind` |
|---|---|---|
| annualised, simple rate | `π`, nominal `rf, rb, re`, real `r, rb_real, re_real` | `:annualised`, `simple = true` |
| level in year h, simple rate | `π`, `rf`, `r` | `:level`, `simple = true` |

Small report helper looping `horizon_percentiles(...; vars = [v])` with title
`"$v, $kind"`. `erp` and `tp` tables are dropped (not in the requested list).

## 10. Fan charts

- Cumulative, wealth of 1 invested (`kind = :wealth`): nominal `rf, rb, re`; real
  `r, rb_real, re_real`.
- Annualised, simple rate (`kind = :annualised, simple = true`): `π`, nominal `rf, rb, re`; real
  `r, rb_real, re_real`.
- Levels, simple rate (`kind = :level, simple = true`, existing chart): `π, r, rf, tp`.

One grid per group, unconditional-mean line from `μ = unconditional_mean(c, A)` mapped to
each series (e.g. `re` → `μ_π + μ_r + erp_target`, real → minus `μ_π`).

## Implementation order

1. `data.jl` real returns + `returns_from_components` rename (+ tests).
2. `io.jl` save/load params, add `YAML`, `CSV` to `Project.toml` (+ round-trip test).
3. `ScenariosEvaluation`: `fan`, `horizon_values(simple = true)`, `path_moments(hist)`,
   `scenarios_correlations` (+ tests). Run `julia --project -e 'using Pkg; Pkg.test()'`.
4. EDA: new vars, fix stale `tp` description in the *Load data* cell (still says "term
   spread"; it is now `rb − rf`), replace *Conclusions* with an *Export* section that loops
   over the three specs and writes `specs/<spec>/`.
5. Report: load spec, remove fitting and data loading, new tables and fans.
6. `render_reports.sh`: loop over `specs/*/`; regenerate `reports/`.

## Acceptance

- `issm_report.qmd` contains no `fit_ar1`, `load_shiller_annual`, `issm_bootstrap_matrix`.
- Rendering the 1935–2025 spec gives the same simulated numbers as before the refactor
  (same seed, same `U`, same `c`, `A`, `y0`, up to float round-trip through YAML/CSV),
  apart from the new/renamed tables.
- Tests pass.
