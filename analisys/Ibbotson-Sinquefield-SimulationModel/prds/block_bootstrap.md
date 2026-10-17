### Block-bootstrap scenarios of real returns

New Quarto report [block_bootstrap_report.qmd](../block_bootstrap_report.qmd) that simulates
real returns by block-bootstrapping historical years, after shifting their means to capital
market assumptions (CMA) given in the report's Quarto params. Evaluation mirrors
[issm_report.qmd](../issm_report.qmd).

## Workflow

```
data/shiller_annual.csv ──▶ block_bootstrap_report.qmd ──▶ reports/bb_<period>_b<block>.html
                             params: cma means, period, block_size, …
                             subperiod → match_means → block_bootstrap → evaluation
```

- Variables: real returns `r` (T-bill), `rb_real`, `re_real`, decimal log units, as in
  `shiller_annual.csv` (written by `export_shiller_annual.jl`). Nominal returns and π are
  out of scope (bootstrapping π alongside would be a later extension).
- Whole historical years are resampled, so the cross-sectional correlation of a year is
  kept exactly; blocks keep autocorrelation up to the block length.
- Units as in the ISSM report: everything in decimal logs, ×100 only in tables and plots;
  annualised values and levels shown as simple rates (`simple = true`).

## 1. Fix and finish `block_bootstrap.jl`

The current file is not included in the module and has several bugs. Target API:

```julia
block_bootstrap(X::AbstractMatrix, n_steps, n_scen; block_size = 5, random_block = false,
                rng = Random.default_rng()) -> Array{Float64,3}          # n_vars × n_steps × n_scen
block_bootstrap(ta::TimeArray, n_steps, n_scen; start_year, kwargs...) -> Scenarios
```

Changes vs. the current draft:

- **Output length `n_steps`** independent of the sample length (the report simulates `T = 100`
  years from ~90 historical ones); today it is forced to the sample length.
- **`rng` keyword** for reproducible reports (today it uses the global RNG).
- **Circular blocks**: start uniformly in `1:n_obs` and wrap indices with `mod1`. The
  current `start ∈ 1:(T − b + 1)` plus truncation at `T` under-samples the first and last
  years and, with random lengths, cuts blocks short. With wrapping every year has the same
  probability, so the bootstrap mean equals the sample mean in expectation — required for
  the CMA match to carry over to the scenarios.
- **Random block lengths** (stationary bootstrap, Politis–Romano): `1 + rand(rng, Geometric(1/block_size))`.
  `Distributions.Geometric` counts failures (support `0, 1, …`), so the current
  `rand(Geometric(1/b))` gives zero-length blocks and mean length `b − 1`. Add `Distributions`
  to the module imports (already a package dependency).
- `block_size ≥ 1` and `block_size ≤ n_obs` checked (`ArgumentError`).
- Fix the docstring (name `block_bootstrap_mv`, wrong `(T, N, B)` output shape).
- `TimeArray` method: uses `colnames(ta)` as `Scenarios.names`; `start_year` kwarg (default
  `year(last(timestamp(ta))) + 1`).
- `overlapping_windows`: not used by the report; either delete it or give it a docstring
  and test. Recommendation: delete (out of scope).
- `include("block_bootstrap.jl")` in `IbbotsonSinquefield.jl`, export `block_bootstrap`.

Tests (`test/ibbotson_sinquefield_test.jl`):

- shape `(n_vars, n_steps, n_scen)` and `Scenarios` names/start_year;
- every simulated year (column) equals some row of `X` (rows are kept whole);
- `block_size = n_obs`, `random_block = false`: each path is a rotation of the sample;
- with fixed blocks, consecutive years inside a block are consecutive sample rows (mod `n_obs`);
- same `rng` seed → identical output;
- mean over many paths ≈ sample mean (circularity), `random_block` mean block length ≈ `block_size`;
- `ArgumentError` for `block_size = 0` and `block_size > n_obs`.

## 2. Report `block_bootstrap_report.qmd`

Quarto params (same `@isdefined` defaults cell as `issm_report.qmd`):

```yaml
params:
  # Capital market assumptions: target means of annual log real returns (decimals)
  cma: {r: 0.005, rb_real: 0.015, re_real: 0.045}
  period: [1935, 2025]      # historical sample bootstrapped
  block_size: 5
  random_block: true        # stationary bootstrap
  T: 100
  horizons: [1, 5, 10, 20, 30, 50, 100]
  n_scen: 10000
  seed: 1234
```

- `cma` holds **means of annual log returns** — the quantity `match_means` shifts. A simple
  expected return `E[R]` maps to roughly `log(1 + E[R]) − σ²/2`; convert before entering it.
- In the defaults cell, `@isdefined(cma) || (cma = Dict("r" => 0.005, …))`; Quarto passes a
  mapping as a `Dict{String}`, so convert once with
  `targets = (; (Symbol(k) => Float64(v) for (k, v) in cma)...)` for `match_means`. Unknown
  names are rejected by `match_means`.
- Other assumption sets: `quarto render … -P cma.re_real:0.05` does not work for nested
  keys, so use `--execute-params <file>.yml` with a full `cma:` mapping, or edit the
  front matter. No separate CMA file or loader in the package.

Sections:

1. **Data and assumptions** — load `shiller_annual.csv`, keep `[:r, :rb_real, :re_real]`,
   `subperiod(..., period)`. Table per variable: historical mean, CMA target,
   shift, all in %. Parameters table (period, block size, random/fixed, T, n_scen, seed).
2. **Historical moments** — `describe_series` and `correlation_table` of the *matched*
   sample (`hist`); same std, skew, autocorrelation and correlations as raw data, means =
   targets. Computed in the report (no spec directory: there is nothing estimated to export).
3. **Simulation** — `sc = block_bootstrap(match_means(hist, targets), T, n_scen; block_size,
   random_block, rng = MersenneTwister(seed))`. Print the achieved mean over all paths vs.
   target. No post-simulation `match_means` on `sc` (would hide bootstrap noise; the circular
   design already makes it unbiased) — open point, see below.
4. **Fan charts** (reuse `fan`, `μ = target`):
   - wealth of 1 invested, `kind = :wealth`: `r, rb_real, re_real`;
   - annualised simple rate, `kind = :annualised, simple = true, scale = 100`.
   Yearly-level fans are dropped (i.i.d.-ish yearly returns, not informative).
5. **Per-path moments** — `path_moments(pct(sc); hist = pct(hist moments))` per variable.
   Note for the reader: within-path `ar1` is biased towards 0 relative to history when
   `block_size` is small (autocorrelation beyond the block is lost).
6. **Correlations** — `mean_path_correlations(sc)` and `path_correlations(sc; hist = hist_cor)`.
7. **Horizon percentiles** — annualised simple rates for `r, rb_real, re_real` at
   `horizons ∩ 1:T` (same `horizon_tables` helper).
8. **Drawdowns** — `drawdown_table` over the first 50 years, max drawdown in %, length in
   years (same `dd_table` helper).

Shared helpers (`show_table`, `pct`, `horizon_tables`, `dd_table`) are copied from
`issm_report.qmd` for now; moving them into `ScenariosEvaluation` is a separate refactor.

`fan` labels the reference line "uncond. mean"; add a `μ_label = "uncond. mean"` kwarg and
pass `"CMA target"` here.

## 3. Rendering

- `render_bb_reports.sh`: render with the front-matter params (extra quarto args passed
  through, e.g. `--execute-params alt.yml`), output `reports/bb_<p1>-<p2>_b<block_size>.html`
  (move file, as in `render_reports.sh`).
- `.gitignore`: ad-hoc `block_bootstrap_report.html` / `_files/`.

## Implementation order

1. `block_bootstrap.jl` fixes + include/export + tests.
2. `fan(...; μ_label)` + smoke test.
3. Report, then render script. Run `julia --project -e 'using Pkg; Pkg.test()'` and render once.

## Acceptance

- Tests pass; `block_bootstrap` is reproducible for a given seed.
- Rendered report: simulated mean of each variable over all years and paths within
  ~2 standard errors of the CMA target; std, skewness and contemporaneous correlations close
  to the matched historical sample.

## Open points

- Match means before bootstrapping only (recommended), or also re-match `sc` exactly?
- Block size: fixed 5, or a short sensitivity table (1, 5, 10) in the report?
- Real-only, or bootstrap `π` with the returns to also produce nominal returns and inflation?
