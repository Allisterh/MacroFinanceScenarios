### Δlog CAPE — EDA

Notebook `analisys/delta_cape_eda.ipynb`: exploratory analysis of the annual valuation change `dlcape = Δlog CAPE` and its relation to the ISSM macro building blocks.

**Goal:** decide how `dlcape` should be treated in the scenarios: as i.i.d. noise around zero (the current ISSM assumption, where it is stripped from erp), as a function of the macro state, or as mean-reverting through the CAPE level.

#### Data and conventions (decisions made)

- Source: `data/shiller_annual.csv` (written by `export_shiller_annual.jl`), annual Dec→Dec, 1935–2025, log units. Load and scale exactly as in `issm_eda.ipynb`: everything × 100 (percent) except the `cape` level.
- Variables used:
  - `dlcape`: target
  - `π`: inflation
  - `r`: ex-post real short rate `rf − π`
  - `tp`: term spread (10y yield − T-bill yield) at the **start** of the year, so it is known at t−1
  - `fr`: fundamental return `re − dlcape`
  - `logcape_lag = log CAPE_{t−1}`: added. The level is the standard valuation predictor, and differences alone cannot show mean reversion.
- Setup cell copied from `issm_eda.ipynb` (`Pkg.activate`, `show_table`). Reuse `describe_series`, `correlation_table`, `subperiod` and `fit_ar1`.
- No new package dependencies, in line with `ar1.jl`: OLS via `X \ y`, and ADF, KPSS and Newey–West written as small helpers in the notebook, like `ljung_box` in `issm_eda.ipynb`.
- Sample periods: the same as in `issm_eda.ipynb` (`full` 1935–2025, `pre1975`, `post1975`), so results line up with the ISSM specs.

#### Caveats to state in the notebook (checked on the data)

- **`fr` is mechanically tied to `dlcape`.** `fr ≡ re − dlcape` holds exactly. Roughly, `fr ≈ dividend yield + π + Δlog(10y-avg real earnings)`, while `dlcape = Δlog real price − Δlog(10y-avg real earnings)`. The earnings-growth term enters both with opposite signs, so part of corr(dlcape, fr) = −0.42 is built in, not economic.
- **`fr` is nominal:** corr(fr, π) = 0.71. Also show the real version `fr_real = fr − π` (corr with dlcape is −0.23).
- **`r` contains `−π`**, so `π` and `r` are mechanically negatively correlated. Expect collinearity in the joint regression and do not read their coefficients in isolation.
- **Contemporaneous vs. predictive.** Same-year regressors (`π`, `r`, `fr`) *explain* dlcape but are not known in advance. `tp`, `logcape_lag` and `dlcape_lag` are predictive. For the simulation the distinction matters: contemporaneous links can be kept through the joint bootstrap, and predictive ones need a model term.

#### Tasks

1. **Load.** Build the TimeArray with `dlcape`, `π`, `r`, `tp`, `fr`, `fr_real`, `logcape` and their lags. Plot `dlcape` alongside `log CAPE` (the cumulated `dlcape`).
2. **Moments.** Show the `describe_series` table of `dlcape` per period (mean, std, skew, excess kurtosis, AR(1), % negative, quantiles). Add a column with the t-statistic of the mean `mean / (std/√n)`: is the historical re-rating significantly different from 0? List the five largest absolute years (e.g. 1937, 1974, 2008).
3. **Correlations.** Correlate `dlcape` with `π`, `r`, `tp`, `fr`, `fr_real` in two blocks, contemporaneous and with the regressors lagged one year, per period (`correlation_table`). Mark |ρ| > 1.96/√n as significant.
4. **Moving averages.** Plot `dlcape` with **trailing** 10y and 20y moving averages (the first 9/19 years are missing) and a zero line. Optional: a 20y rolling std, to check whether the volatility is stable.
5. **Stationarity.** Test both `dlcape` and `log CAPE`:
   - ADF with a constant. Choose the lag by BIC up to 4 lags and use the MacKinnon 5% critical value for n≈90 (≈ −2.89).
   - KPSS (level-stationary null, Newey–West bandwidth) as the complementary test.
   - Report one table: statistic, critical value and decision per series and period.
   - Interpretation: dlcape should be clearly stationary. Whether log CAPE is stationary (H0 unit root is likely not rejected in 90 years) decides if a mean-reversion term `logcape_lag` is defensible.
6. **Univariate regressions.** Fit `dlcape_t = a + b x_t + e_t`, one regression for each of `π`, `r`, `tp`, `fr`, `fr_real`, `logcape_lag`. Show one comparison table: b, Newey–West t-stat (lag 2), R², adjusted R², n, per period.
7. **Joint regression and tests of impact.**
   - Model A (contemporaneous): `dlcape_t = a + β_π π_t + β_r r_t + β_tp tp_t + β_fr fr_t + ρ dlcape_{t−1} + e_t`
   - Model B (predictive only): `dlcape_t = a + β_tp tp_t + γ logcape_{t−1} + ρ dlcape_{t−1} + e_t`
   - Report R² and adjusted R², with residual diagnostics: `ljung_box` and a residual plot.
   - Tests of impact:
     - **Individual:** t-tests with Newey–West (HAC) standard errors. Residuals are likely heteroskedastic (large crisis years).
     - **Joint:** a Wald/F-test that all macro coefficients are 0 (H0: dlcape is unpredictable, i.e. the i.i.d. assumption).
     - **Incremental:** drop one variable at a time and report ΔR² with a nested F-test (partial F). Because π, r and fr are collinear, this is more informative than the individual t-tests.
     - **Stability:** a Chow test at 1975, with the `pre1975` and `post1975` coefficients side by side.
     - **Out-of-sample** (for Model B, the only one usable for forecasting): expanding-window R²_OS against the historical mean (Campbell–Thompson), starting in 1965.
8. **Conclusion cell.** Summarise which of the three treatments of `dlcape` the evidence supports, and what that means for the ISSM simulation.

#### Out of scope

- No changes to `src/` or to the ISSM specs or report. If a model for dlcape is adopted, that is a separate plan.
