# Trend-Cycle Macro-Finance Model: Asset-Pricing Justification

Oct 27, 2026 · @Mat

## Verdict

Yes, the model is appropriate, provided six specification details are handled correctly. Writing the dividend yield as short-rate trend + equity-premium trend − dividend-growth trend + cycle is the Gordon growth model applied to stochastic trends, or equivalently the Campbell–Shiller (1988) present-value identity in which random-walk trends are their own long-run forecasts. It is a reduced-form, identity-consistent model, not a structural equilibrium model, and it sits in an established literature: present-value state-space models estimated with the Kalman filter (van Binsbergen and Koijen 2010), shifting steady-state valuation ratios (Lettau and Van Nieuwerburgh 2008), and trend-cycle VARs for interest rates (Del Negro, Giannone, Giannoni and Tambalotti 2017).

What the decomposition buys you: in the long run the expected equity return equals short-rate trend + ERP trend, with dividend growth cancelling exactly, and a shock to any trend produces the correct immediate revaluation of prices.

What must be right for it to hold (detailed in *Weaknesses and recommended changes*):

1. Short rate and dividend growth in the same units: both nominal, or both real with an explicit inflation trend.
2. Total payout (dividends + net buybacks), or a payout trend; otherwise the ERP trend absorbs the post-1980 shift to buybacks.
3. Enter the yield as log(1 + D/P): unit loadings are then exact, and the model stays linear.
4. Returns computed from the price path (yield + growth + revaluation), never as yield + growth alone.
5. Tight priors on trend-shock variances, because random-walk trends make 30-year simulations explode and the ERP trend is identified only as a residual.
6. Treat the random-walk trends as a local approximation: theory requires r + ERP − g > 0 and a stationary valuation ratio.

## Model specification

Each observable is a constant plus a linear combination of random-walk trends plus one stationary cycle; the cycles jointly follow a VAR. All rates are per period, continuously compounded, in decimals (divide annual rates by 4 at quarterly frequency).

```latex
\text{Observation:}\quad y_t = a + \Lambda\,\tau_t + c_t \;(+\,u_t)
```

```latex
\text{Trends:}\quad \tau_t = \tau_{t-1} + \eta_t,\qquad \eta_t \sim N(0, Q_\tau)
```

```latex
\text{Cycles:}\quad c_t = \Phi_1 c_{t-1} + \dots + \Phi_p c_{t-p} + \varepsilon_t,\qquad \varepsilon_t \sim N(0, \Sigma_c)
```

The state vector stacks the trends and p lags of the cycles, so the system is linear and Gaussian and the Kalman filter and smoother give exact likelihood and state estimates. The optional measurement error u\_t is usually dropped; then each cycle equals its observable minus its trend part.

Recommended loading matrix Λ (rows = observables, columns = trends; blank = 0). Inflation and a real-rate trend are added to your four asset blocks so that units stay consistent (see the units section).

| Observable | Inflation trend τπ | Real short-rate trend τr | Term-premium trend τtp | Real payout-growth trend τg | ERP trend τerp | Commodity trend τs | Cycle |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Inflation π | 1 |  |  |  |  |  | cπ |
| Short rate i (3-month bill) | 1 | 1 |  |  |  |  | ci |
| Long yield y (10-year) | 1 | 1 | 1 |  |  |  | cy |
| Nominal payout growth Δd | 1 |  |  | 1 |  |  | cg |
| Payout yield x = log(1 + D/P) | 0 (= 1 − 1) | 1 |  | −1 | 1 |  | cx |
| Commodity: real spot change or roll yield | (1 if nominal) |  |  |  |  | 1 | cs |

The x row is your specification, with the yield entered as x = log(1 + D/P): short-rate trend + ERP trend − nominal growth trend + one cycle. The unit loadings are exact for the Lettau–Van Nieuwerburgh steady state, derived in the next section. Because the nominal short rate and nominal growth both load on τπ, inflation cancels out of the valuation ratio.

## Asset-pricing justification of the equity block

The specification follows in four steps from an accounting identity plus one assumption (no rational bubbles); no particular utility function or equilibrium is needed.

**Step 1 — Return identity.** The gross return is R = (P' + D')/P. Log-linearising around the mean log price–dividend ratio (Campbell and Shiller 1988):

```latex
r_{t+1} \approx k + dp_t - \rho\, dp_{t+1} + \Delta d_{t+1},\qquad \rho = \frac{1}{1+\overline{D/P}} \approx 0.96\text{--}0.97 \text{ (annual)}
```

So the return is dividend growth + yield + a revaluation term. This is why equities split naturally into a growth block and a yield block.

**Step 2 — Present value.** Solving forward and ruling out bubbles (ρ^j dp\_{t+j} → 0):

```latex
dp_t = -\frac{k}{1-\rho} + \sum_{j\ge 0} \rho^{j}\, E_t\big[r_{t+1+j} - \Delta d_{t+1+j}\big]
```

The yield is high when expected returns are high or expected growth is low; nothing else can move it.

**Step 3 — Split expected return and growth into trend and cycle.** Write the expected return as short rate + equity premium, E r = i + erp, and each component as trend + cycle. A random walk forecasts itself, E\_t τ\_{t+j} = τ\_t for every horizon, so the trend part of the infinite sum is a geometric series:

```latex
dp_t = a + \frac{1}{1-\rho}\Big(\tau^{i}_t + \tau^{erp}_t - \tau^{g}_t\Big) + c^{dp}_t,\qquad c^{dp}_t = \sum_{j\ge0}\rho^j E_t\big[c^{i} + c^{erp} - c^{g}\big]_{t+j}
```

This is your specification for the log yield, with loading 1/(1−ρ) ≈ 25–30 and a = −k/(1−ρ). It is accurate only near the mean yield; Step 4 gives an exact linear alternative.

**Step 4 — An exact linear form: observe log(1 + D/P).** Lettau and Van Nieuwerburgh (2008) write the steady state exactly, as Gordon growth in gross terms: D/P = (e^μ − e^g)/e^g, where μ = τi + τerp is the log expected gross return and g = τg. Rearranging gives a relation that is linear in the trends with unit loadings:

```latex
x_t \equiv \log\!\Big(1 + \frac{D_t}{P_t}\Big) = \tau^{i}_t + \tau^{erp}_t - \tau^{g}_t + c^{x}_t
```

The nonlinearity lies only in the map from trends to the observable, so transforming the data before filtering keeps the Kalman filter linear and exact. The identity holds for the trailing yield D\_t/P\_t, the series Shiller and index providers publish. Prices are recovered as P = D/(e^x − 1), which is positive exactly when short rate + ERP > growth, the transversality condition.

The transform removes the approximation error of the alternatives. Price error versus the exact Lettau–Van Nieuwerburgh mapping, linearised at a 3% yield:

| True D/P | Plain D/P, unit loadings | Log D/P, loadings 1/(1−ρ) ≈ 34 | log(1 + D/P), unit loadings |
| --- | --- | --- | --- |
| 1% | +0.5% | −35% | 0 |
| 2% | +1.0% | −7% | 0 |
| 4% | +2.0% | −4.5% | 0 |
| 6% | +3.1% | −26% | 0 |

"Exact" refers to the steady-state mapping, as in Lettau and Van Nieuwerburgh. With stochastic trends, small convexity terms from future trend shocks remain; they are second order when trend variances are small and are absorbed in the constant.

**Two consequences that make the model coherent for simulation:**

- Long-run expected return = τi + τerp exactly. Prices grow at e^g, so the gross return is e^g · (1 + D/P) = e^g · e^(μ−g) = e^μ. Dividend growth cancels: faster growth is offset one-for-one by a lower yield, so a growth assumption cannot inflate long-run returns.
- Trend shocks reprice equities immediately. A permanent +1pp shock to the short-rate trend raises x from 3% to 4%, so P/D falls from 32.8 to 24.5: prices drop about 25% at once, then expected returns are 1pp higher forever. P/D is equity duration, about 30 years.

Because μ is the log of the expected gross return, this ERP is close to the arithmetic premium and already includes σ²/2. The expected log excess return is lower by roughly 1.2–1.5pp at 15–17% volatility; report which one you quote.

## Discount anchor and units

The short-rate trend is the correct anchor, and inflation must enter the short rate and dividend growth identically.

**Why the short rate, not the bond yield.** The present-value sum discounts at the whole future path of expected one-period returns. Under a random walk, the trend of the short rate is the market's expected average future short rate at every horizon, which is exactly the infinite-horizon riskless discount rate. Using the 10-year yield instead would add the bond term premium to the equity discount rate, double-counting duration risk and making the ERP a premium over bonds rather than over cash. Keeping the ERP over cash also makes stocks and bonds comparable: both earn the short-rate trend plus their own premium trend.

**Nominal vs real.** If the short rate is nominal, dividend growth must be nominal too. Then the inflation trend enters both with loading +1 and cancels in the yield, as in the loading table: long-run equity valuations are neutral to inflation, as present-value logic requires. Mixing a nominal short rate with real growth reproduces the "Fed model" error, comparing a nominal yield to a real yield, which Modigliani and Cohn (1979) called money illusion and Asness (2003) criticised. Decomposing the short-rate trend into r\* + π\* has two further payoffs: the real-rate trend can be anchored to r\* estimates, and simulated real returns for all four assets come out directly.

**Short-run inflation effects stay in the cycles.** Empirically, stocks fall when inflation rises unexpectedly (Fama 1981). The model allows this through the VAR on the cycles without breaking long-run neutrality.

## The single dividend-yield cycle

One unsplit cycle for the yield is consistent with theory and loses nothing, because the split is recovered from the VAR after estimation.

- **The yield is a sufficient statistic.** Step 2 shows the yield cycle is the present value of cyclical expected returns minus cyclical expected growth. Prices reveal only that difference, so splitting the cycle a priori would add labels the data cannot separate.
- **The split falls out of the VAR.** With returns built from the identity, the one-step expected return is E\_t r = k + dp\_t − ρE\_t dp\_{t+1} + E\_t Δd\_{t+1}, all VAR forecasts. Discount-rate and cash-flow news then follow from Campbell (1991). Expect the yield cycle mainly to forecast returns, not dividend growth (Cochrane 2008).
- **No inconsistency can arise** as long as there is no separate equity-return equation. Equity returns are an output of the yield and growth blocks, never an input.
- **Optional structural upgrade.** Make the expected-return cycle and expected-growth cycle latent AR(1) states and let the yield cycle be their implied present value (van Binsbergen and Koijen 2010). This stays linear and Kalman-filtered, uses fewer parameters, and its cross-equation restrictions can be tested against your unrestricted VAR by likelihood ratio.

## Bonds and commodities in the same framework

Bonds follow the same present-value logic as equities over a finite horizon; commodities have no cash flows, so their block rests on futures-return accounting instead.

**Bonds.** The expectations hypothesis plus a term premium gives the long yield as the average expected short rate plus the premium. A random-walk short-rate trend is its own forecast at every maturity, so:

```latex
y^{(n)}_t = \tau^{i}_t + \tau^{tp}_t + c^{y}_t,\qquad r^{bond}_{t+1} \approx y^{(n)}_t - D_n\big(y^{(n-1)}_{t+1} - y^{(n)}_t\big) + \tfrac12 C_n \big(\Delta y\big)^2
```

Here D\_n and C\_n are duration and convexity of the constant-maturity bond, and y^(n−1) should include roll-down if you model the curve slope. The long-run expected bond return is short-rate trend + term-premium trend, mirroring equities.

The shared trends give a stock–bond correlation that depends on the source of the shock. A real-rate trend shock hits both assets through duration: positive correlation of returns. An inflation trend shock hits bonds but cancels out of equity valuations: low correlation. The cycle covariance adds the short-run effects, so the model can produce both the positive correlation of the 1970s–1990s and the negative one of 2000–2020.

**Commodities.** A fully collateralised futures position earns:

```latex
r^{cmdty}_{t+1} = i_t + \underbrace{\Delta s_{t+1}}_{\text{spot change}} + \underbrace{roll_{t+1}}_{\text{basis / roll yield}}
```

The excess return over cash is therefore spot change + roll yield, and the theory of storage (Kaldor, Working) and hedging pressure (Keynes) explains the roll yield (Gorton and Rouwenhorst 2006; Erb and Harvey 2006). Two specification rules follow:

- Put the trend on the roll yield and an inflation loading on the nominal spot change. Avoid a random-walk trend on the growth rate of real spot prices: it makes log real prices I(2), whereas over very long samples real commodity prices show little or no drift.
- Collateral return is the short rate, so commodities automatically share the short-rate trend with the other assets.

## Weaknesses and recommended changes

The most important weakness is that the ERP trend is identified only as a residual, so it absorbs every other misspecification; the table runs from most to least serious.

| Issue | Why it matters | Recommended fix |
| --- | --- | --- |
| ERP trend is a residual | With one trend per observable, τerp = τdp/b − τi + τg. The Gordon restriction is a labelling, not a testable constraint, and any error in the yield or growth data lands in the ERP. | Tight prior on its shock variance. Add observables that load on the same trends: earnings yield or CAPE, survey return expectations, analyst long-term growth. Test restricted versions by likelihood, e.g. a constant ERP implies dp − b(i − g) is stationary. |
| Buybacks | US dividend yields fell structurally after about 1980 as payouts shifted to repurchases (Fama and French 2001). The model reads this as a falling ERP. | Use net payout yield and net payout growth (Boudoukh et al. 2007), or add a payout-ratio trend to the yield row. |
| Random-walk trends at 30+ years | Trend variance grows linearly with horizon. Simulated rates drift far negative or high, and r + ERP − g can turn negative, so x ≤ 0 and P = D/(e^x − 1) is undefined. Theory needs a stationary valuation ratio. | Small trend variances via priors, or trends with roots like 0.99 for simulation, so that x stays positive. Report the share of paths breaching r + ERP > g. Do not clip paths, which breaks the identities. |
| Pile-up of trend variances | Maximum likelihood often estimates random-walk variances at zero (Stock and Watson 1998). | Bayesian estimation with informative priors on Qτ (as in Del Negro et al. 2017), or median-unbiased estimation. |
| Trend vs cycle in the yield | The yield cycle is very persistent, so a slow cycle and a trend shift look alike (Lettau and Van Nieuwerburgh 2008). | Long samples (Shiller data from 1871; the Kalman filter handles missing series), a prior capping VAR persistence, and sensitivity tables. |
| I(2) dividends | A random-walk trend in growth makes log dividends I(2). Long-run payout growth is tied to GDP growth less net dilution (Bernstein and Arnott 2003). | Very tight prior on τg, or load it on a GDP-growth trend observable. |
| Linearisation error | Using log(1 + D/P) removes linearisation error in the yield. The Campbell–Shiller identity still used for return decompositions degrades far from the mean yield, e.g. a 1.1% yield in 2000. | Compute simulated returns from prices, not the identity; for decompositions, set ρ from the sample mean and check the identity residual. |
| Gaussian shocks | Crash frequency and volatility clustering are understated, so long-horizon tails are too thin. | Keep Kalman estimation; simulate with bootstrapped or Student-t residuals. |
| Zero lower bound | A linear Gaussian model cannot keep rates non-negative. | Accept it, or use a shadow short rate (e.g. Wu–Xia) as the observable in ZLB periods. |

## Estimation and simulation

Estimate the state-space model by Kalman-filter likelihood with informative priors on trend variances, then simulate levels and prices and compute returns exactly from them.

**Estimation**

1. Data, quarterly: 3-month bill, 10-year government yield, CPI inflation, index dividends or net payout and the payout yield (Shiller's series give a history from 1871), and a commodity index's spot and excess returns. The Kalman filter treats series that start later as missing values.
2. Initial states: diffuse or a wide prior for the trends; the stationary distribution for the cycles.
3. Likelihood from the prediction-error decomposition. Use either maximum likelihood with fixed signal-to-noise ratios for the trends, or Bayesian estimation: inverse-Wishart prior with a small scale on Qτ, a Minnesota prior on the VAR, and Gibbs sampling with the Durbin–Koopman simulation smoother. Both keep inference within the linear Kalman framework.
4. Validation: compare smoothed trends with external estimates (Holston–Laubach–Williams r\*, ACM term premium, Damodaran's implied ERP); check the return-identity residual; test out-of-sample return forecasts against the historical mean.

**Simulation**

1. Draw parameters from the posterior (or fix them) and the starting state from the smoothed distribution at the last date, N(τ\_T|T, P\_T|T). Trend uncertainty at the start dominates 30-year dispersion, so do not start from a point estimate.
2. Simulate trends and cycles forward and rebuild every observable from the loading matrix.
3. Build levels, then returns:
   - Equities: log D accumulates Δd; P = D/(e^x − 1); total return = log(P' + D') − log P. This includes the revaluation term automatically.
   - Bonds: the duration–convexity formula above on the simulated yield path.
   - Cash: the simulated short rate.
   - Commodities: short rate + spot change + roll yield.
4. Report annualised return quantiles by horizon, cross-asset correlations, and the share of paths breaching r + ERP > g or the zero bound.

Because x enters the measurement equation linearly, the steady-state expected log return is exactly short-rate trend + ERP trend at every date and horizon. The full distribution of cumulative returns needs simulation, since P = D/(e^x − 1) is nonlinear; use the steady-state identity to check the simulation.

## References

Cited from memory, not looked up for this draft; verify volume and page details before external use.

- Adrian, Crump and Moench (2013). Pricing the term structure with linear regressions. *Journal of Financial Economics* 110.
- Asness (2003). Fight the Fed model. *Journal of Portfolio Management* 30(1).
- Bernstein and Arnott (2003). Earnings growth: the two percent dilution. *Financial Analysts Journal* 59(5).
- Boudoukh, Michaely, Richardson and Roberts (2007). On the importance of measuring payout yield. *Journal of Finance* 62.
- Campbell (1991). A variance decomposition for stock returns. *Economic Journal* 101.
- Campbell and Shiller (1988). The dividend-price ratio and expectations of future dividends and discount factors. *Review of Financial Studies* 1(3).
- Cochrane (2008). The dog that did not bark: a defense of return predictability. *Review of Financial Studies* 21(4).
- Del Negro, Giannone, Giannoni and Tambalotti (2017). Safety, liquidity, and the natural rate of interest. *Brookings Papers on Economic Activity*, Spring.
- Durbin and Koopman (2002). A simple and efficient simulation smoother for state space time series analysis. *Biometrika* 89.
- Erb and Harvey (2006). The strategic and tactical value of commodity futures. *Financial Analysts Journal* 62(2).
- Fama (1981). Stock returns, real activity, inflation, and money. *American Economic Review* 71.
- Fama and French (2001). Disappearing dividends. *Journal of Financial Economics* 60.
- Gorton and Rouwenhorst (2006). Facts and fantasies about commodity futures. *Financial Analysts Journal* 62(2).
- Holston, Laubach and Williams (2017). Measuring the natural rate of interest: international trends and determinants. *Journal of International Economics* 108.
- Lettau and Van Nieuwerburgh (2008). Reconciling the return predictability evidence. *Review of Financial Studies* 21(4).
- Modigliani and Cohn (1979). Inflation, rational valuation and the market. *Financial Analysts Journal* 35.
- Stock and Watson (1998). Median unbiased estimation of coefficient variance in a time-varying parameter model. *Journal of the American Statistical Association* 93.
- van Binsbergen and Koijen (2010). Predictive regressions: a present-value approach. *Journal of Finance* 65(4).
- Wu and Xia (2016). Measuring the macroeconomic impact of monetary policy at the zero lower bound. *Journal of Money, Credit and Banking* 48.
