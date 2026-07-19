# Build Specification — Empirical Experiment Planner + min-n-per-arm Power Calculator

**Version:** v01, 2026-07-19
**Author:** Flaudisegna (Claude Fable 5) — specification and verification
**Builders:** Claude Sonnet subagents (Builder A: stats engine; Builder B: app shell)
**Plan of record:** `docs/plan-v2-max-effort-improvements_2026-07-19.md` (v2 supersedes v1; v1 preserved for the boundary principle's origin)

This spec is law for the build. Builders MUST NOT silently deviate. If you believe the spec contains an error, STOP work on that item and report the discrepancy in your final report — do not "fix" pinned values or formulas on your own judgment.

---

## 1. Product overview

A single self-contained HTML file: an interactive planner for controlled empirical experiments with an exact min-n-per-arm power calculator. Two halves, one data model:

- **Planner** — research question, hypotheses, IVs/DVs, auto-generated condition grid (factorial detection), randomization plan with a working seeded generator, confound register, validity-threat checklists, analysis plan, exports.
- **Calculator** — closed-form power/sample-size for the full test family in §6, exact noncentral-distribution numerics (§5), derive-from-what-you-know effect-size helpers, MDE inversion, precision mode, cluster and attrition adjustments, live self-validating receipts panel (§7).

### 1.1 Principles (these shape copy and behavior)

1. **The closed-form boundary.** Everything with a closed-form, assumption-light formula is in. Everything requiring an invented data-generating model is out, and the tool SAYS SO and points to simulation (§11.8). No fake precision.
2. **Receipts or it didn't happen.** The engine validates itself against published anchors live in the page, pass/fail visible (§7).
3. **You never pick an abstract effect size.** Every calculator has a derive helper: means→d, proportions→h, cell means→f per term, paired inputs→dz (§6.14).
4. **SESOI doctrine.** Power for the smallest effect of interest, not the hoped-for effect. Pilot > literature > benchmarks; Cohen benchmarks labeled "last resort" (§11.2).
5. **Unit-agnostic.** A unit label ("participants", "prompts", "sessions", "runs") propagates through every generated sentence.
6. **Honest labels everywhere.** Every result card names its test, method, and approximation status (§10.6).

## 2. Delivery constraints (hard)

- ONE self-contained HTML file at repo root: `index.html`. Zero network requests, zero CDN, zero external fonts/images. Opens offline via `file://` double-click.
- Dark background + bright neon accents (tokens §10.2). Body text is bright off-white — neon is for accents, never body copy.
- **Vertigo-safe motion (hard accessibility rule):** NO rotating, spinning, or spiral animation or iconography anywhere — no spinner loaders, no spiral emoji/glyphs, no rotating hourglass. Transitions limited to opacity/translate. Respect `prefers-reduced-motion: reduce` (disable all transitions). This is a non-negotiable user accessibility requirement.
- Colorblind-safe: accents differ in lightness as well as hue; meaning is never encoded by color alone (always text/shape too).
- Keyboard operable; visible `:focus-visible` outlines; inputs labeled with `<label for>`; result regions `aria-live="polite"`; hit targets ≥ 40px.
- Vanilla JS/CSS only. No frameworks, no build-time dependencies beyond node built-ins.
- All randomness in the randomizer flows from the user-visible seed field (deterministic). Never call `Date.now()`/`Math.random()` for anything user-facing except a "suggest random seed" button MAY use `Math.random` to fill the seed field (the generated list still depends only on the field value).

## 3. Repository architecture and assembly

```
app-empirical-experiment-planner/
├── index.html                 # BUILT artifact (assembled; committed)
├── build.mjs                  # assembly: inlines engine into shell → index.html
├── src/
│   ├── stats-engine.js        # Builder A. Pure logic, zero DOM. Dual browser/node export.
│   └── app-shell.html         # Builder B. Full app, engine EXCLUDED, marker included.
├── test/
│   └── receipts-runner.cjs    # Builder A. Node receipt runner, zero deps, exit 1 on fail.
└── docs/                      # plans, this spec, audits
```

- `src/app-shell.html` contains, as its FIRST script tag, exactly:
  `<script>/*__INLINE_STATS_ENGINE__*/</script>`
  followed by separate `<script>` tag(s) with all app code. `build.mjs` (already written by the orchestrator) replaces the marker line with the full engine source to produce `index.html`.
- Engine export pattern (end of `stats-engine.js`):
  ```js
  if (typeof module !== "undefined" && module.exports) { module.exports = StatsEngine; }
  if (typeof globalThis !== "undefined") { globalThis.StatsEngine = StatsEngine; }
  ```
- App code calls only the public API in §4. Builder B must not reimplement any math (exception: trivial display arithmetic like totals and percentages).
- Builders work only on their own files. Neither builder commits.

## 4. Stats engine public API (contract — lock signatures)

`StatsEngine` is a plain object. All functions pure and deterministic. Angles/probabilities in natural units. `opts` fields with defaults: `alpha=0.05`, `tails=2` (1 or 2), `power=0.80`, `ratio=1` (ratio = n2/n1).

```js
// ——— distributions (exposed for receipts/identity checks) ———
normCdf(z); normInv(p); erf(x);
tCdf(t, df); tInv(p, df);
fCdf(x, d1, d2); fInv(p, d1, d2);
chisqCdf(x, df); chisqInv(p, df);
nctCdf(t, df, ncp);                       // noncentral t
ncfCdf(x, d1, d2, lambda);                // noncentral F
ncchisqCdf(x, df, lambda);                // noncentral chi-square

// ——— two-arm & one-arm means (exact noncentral t) ———
powerT2(n1, d, opts) -> {power}                    // independent two-sample
solveNT2(d, opts) -> res                           // res shape in §4.1
mdeT2(n1, opts) -> {d}
powerT1(n, d, opts) -> {power}                     // one-sample AND paired (d = dz)
solveNT1(d, opts) -> res
mdeT1(n, opts) -> {d}

// ——— two-arm proportions (normal approx; method-labeled) ———
// method: 'pooled' (default) | 'unpooled' | 'arcsine' ; cc: continuity correction (pooled/unpooled only; arcsine requires ratio=1)
powerProp2(n1, p1, p2, opts&{method, cc}) -> {power}
solveNProp2(p1, p2, opts&{method, cc}) -> res
mdeProp2(n1, p1, opts&{method, cc, direction:+1|-1}) -> {p2}

// ——— shifted-null designs (normal approx) ———
solveNNIMeans(deltaTrue, margin, sd, opts) -> res        // margin > 0; reject inferiority
powerNIMeans(n1, deltaTrue, margin, sd, opts) -> {power}
solveNTostMeans(deltaTrue, margin, sd, opts) -> res      // requires |deltaTrue| < margin
powerTostMeans(n1, deltaTrue, margin, sd, opts) -> {power}
solveNNIProp(p1, p2, margin, opts) -> res                // unpooled SE
solveNTostProp(p1, p2, margin, opts) -> res

// ——— paired proportions (McNemar, Connor approx) ———
solveNMcNemar(p1, p2, rho, opts) -> res                  // derives p01/p10; validates rho range
solveNMcNemarDirect(p10, p01, opts) -> res
mcnemarRhoBounds(p1, p2) -> {lo, hi}

// ——— k-group omnibus (exact noncentral F / chi-square) ———
powerAnova(nPerGroup, f, k, opts) -> {power}
solveNAnova(f, k, opts) -> res                           // res.nPerGroup, res.Ntotal
solveNChisq(w, df, opts) -> res                          // res.Ntotal
solveNChisqKprops(props, opts) -> res                    // k×2 omnibus; also returns {w}
mdeAnova(nPerGroup, k, opts) -> {f}

// ——— factorial (exact noncentral F per term) ———
factorialEffectsFromCells(grid, sd) -> {fA, fB, fAB, rowEffects, colEffects, interactionEffects, grandMean}
solveNFactorialTerm(fTerm, dfTerm, cells, opts) -> res   // res.nPerCell, res.Ntotal; error df = N - cells
powerFactorialTerm(nPerCell, fTerm, dfTerm, cells, opts) -> {power}

// ——— correlation (Fisher z) ———
solveNCorr(r, opts) -> res ; powerCorr(n, r, opts) -> {power} ; mdeCorr(n, opts) -> {r}

// ——— precision mode (CI half-width; exact-t iteration where t applies) ———
nMeanCI(sd, halfWidth, conf) ; nMeanDiffCI(sd, halfWidth, conf)          // per arm
nPropCI(p, halfWidth, conf) ; nPropDiffCI(p1, p2, halfWidth, conf)       // per arm

// ——— adjustments & derivation helpers ———
deff(m, icc, cvm=0) -> number                            // 1 + ((1+cvm^2)*m - 1)*icc
applyCluster(n, m, icc, cvm=0) -> {nAdjusted, clustersPerArm, deff, fewClustersWarning}  // warn if clustersPerArm < 8
applyAttrition(n, rate) -> nInflated                     // ceil(n / (1 - rate))
comparisonsCount(k, scheme) -> m                         // 'allPairs': k(k-1)/2 ; 'vsControl': k-1
bonferroniAlpha(alpha, m) -> alpha/m
rankInflate(n) -> ceil(n / 0.955)                        // ARE vs t under shift alternatives
dFromMeans(m1, m2, sd) ; dFromMeansSds(m1, m2, sd1, sd2) // pooled sd = sqrt((sd1^2+sd2^2)/2)
hFromProps(p1, p2)                                       // 2*(asin(sqrt(p1)) - asin(sqrt(p2))), abs
dzFromPaired(delta, sd, rho)                             // dz = delta / (sd*sqrt(2*(1-rho)))
anovaFFromMeans(means, sd)                               // sqrt(mean((mu_i - mean)^2)) / sd

// ——— receipts ———
runReceipts() -> [{id, desc, expected, got, pass, method, source}]
```

### 4.1 `res` shape for all `solveN*`

```js
{ nPerArm (or nPerGroup/nPerCell/Ntotal as natural), Ntotal, achievedPower, // at returned integer n
  df, criticalValue, ncpAtSolution, method, notes[] }
```
Integer n via: find smallest integer with power ≥ target (double-then-binary-search; power is monotone in n). Per-arm minimum 2 (one-sample/paired minimum 4 with a "tiny-sample" note). If n would exceed 10^8, throw a descriptive error.

## 5. Numerics layer (algorithms — implement exactly)

Target: absolute CDF error ≤ 1e-9 in central distributions, ≤ 1e-8 in noncentral. All receipts in §7 must pass with these algorithms.

- **`erf(x)`** — Abramowitz & Stegun 7.1.26: `t = 1/(1+0.3275911*x)`; `erf ≈ 1 − (a1*t + a2*t² + a3*t³ + a4*t⁴ + a5*t⁵)·exp(−x²)` for x≥0, odd extension for x<0. Coefficients: a1=0.254829592, a2=−0.284496736, a3=1.421413741, a4=−1.453152027, a5=1.061405429. `normCdf(z) = 0.5*(1+erf(z/√2))`.
- **`normInv(p)`** — Acklam's algorithm, coefficients exactly:
  - a = [−3.969683028665376e+01, 2.209460984245205e+02, −2.759285104469687e+02, 1.383577518672690e+02, −3.066479806614716e+01, 2.506628277459239e+00]
  - b = [−5.447609879822406e+01, 1.615858368580409e+02, −1.556989798598866e+02, 6.680131188771972e+01, −1.328068155288572e+01]
  - c = [−7.784894002430293e-03, −3.223964580411365e-01, −2.400758277161838e+00, −2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
  - d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
  - break points p_low=0.02425, p_high=1−0.02425; rational forms per published algorithm.
- **`lgamma(x)`** — Lanczos g=7, n=9: coefficients 0.99999999999980993, 676.5203681218851, −1259.1392167224028, 771.32342877765313, −176.61502916214059, 12.507343278686905, −0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7; reflection formula for x<0.5.
- **Regularized incomplete beta `ibeta(x,a,b)`** — continued fraction (Lentz), symmetry switch when `x > (a+1)/(a+b+2)`: `ibeta = 1 − ibeta(1−x, b, a)`. Front factor `exp(a*ln x + b*ln(1−x) + lgamma(a+b) − lgamma(a) − lgamma(b)) / a`. CF terms: even `d_{2m} = m(b−m)x/((a+2m−1)(a+2m))`, odd `d_{2m+1} = −(a+m)(a+b+m)x/((a+2m)(a+2m+1))`. Tolerance 1e-14, max 400 iterations.
- **Regularized lower incomplete gamma `pgamma(a,x)`** — series when `x < a+1` (term_{n+1} = term_n * x/(a+n+1), from 1/a), else 1 − continued fraction (Lentz). Tolerance 1e-14.
- **Central CDFs** — `tCdf`: for t≥0, `1 − 0.5*ibeta(df/(df+t²), df/2, 0.5)`; mirror for t<0. `fCdf(x) = ibeta(d1x/(d1x+d2), d1/2, d2/2)`. `chisqCdf(x) = pgamma(df/2, x/2)`.
- **Quantiles** `tInv, fInv, chisqInv` — bracket-doubling + bisection on the CDF to 1e-12 relative, 200 iterations cap. (Simple and robust beats clever here.)
- **`nctCdf(t, df, δ)`** — integrate over the chi-square denominator: `P(T'≤t) = ∫ Φ(t·√(v/df) − δ)·f_{χ²_df}(v) dv`, v from `chisqInv(1e-13, df)` to `chisqInv(1−1e-13, df)`, composite Simpson, refine by doubling intervals until |ΔP| < 1e-10 (start 256, cap 32768 intervals). `f_{χ²}` via `exp((df/2−1)ln v − v/2 − lgamma(df/2) − (df/2)ln 2)`. For df > 1e5 use the normal approximation `Φ((t(1−1/(4df)) − δ)/√(1+t²/(2df)))`.
- **`ncchisqCdf(x, df, λ)`** — Poisson mixture: `Σ_j w_j · pgamma((df+2j)/2, x/2)`, `w_j = e^{−λ/2}(λ/2)^j/j!`. Sum two-sided from `j0 = floor(λ/2)` with log-space weight recurrence; stop each direction when `w_j < 1e-16`. λ=0 → central.
- **`ncfCdf(x, d1, d2, λ)`** — Poisson mixture with constant beta argument: `Σ_j w_j · ibeta(v, (d1+2j)/2, d2/2)`, `v = d1x/(d1x+d2)`, same weight scheme.
- **Power via noncentral distributions:**
  - Two-sample t: df = n1+n2−2, ncp `δ = d·√(n1·n2/(n1+n2))`. Two-tailed: `P(T'>t_{1−α/2}) + P(T'<−t_{1−α/2})` (include the second term). One-tailed: `P(T'>t_{1−α})`.
  - One-sample/paired: df = n−1, `δ = d·√n`.
  - ANOVA/factorial term: `λ = f²·N_total`, power = `1 − ncfCdf(fInv(1−α, df1, df2), df1, df2, λ)`; ANOVA df1=k−1, df2=N−k; factorial term df1=dfTerm, df2=N−cells.
  - χ²: `λ = w²·N`, power = `1 − ncchisqCdf(chisqInv(1−α, df), df, λ)`.

## 6. Calculator formulas (per test)

Shared: `zA = normInv(1−alpha/tails)`, `zP = normInv(power)`. "Solve n" = integer search on the EXACT power function where one exists (t/F/χ² families); normal-approx families solve their formula then verify/adjust by evaluating their stated power function at integer n.

1. **Two-sample t** — exact noncentral t (§5). Allocation ratio r: n2 = ceil(r·n1). Assumes equal variances for planning (copy: §11.9).
2. **One-sample / paired t** — same machinery, df=n−1. Paired UI derives `dz` via `dzFromPaired` (σ_d = σ√(2(1−ρ))).
3. **Two proportions, pooled z (default)** — with r = n2/n1, p̄=(p1+r·p2)/(1+r), q̄=1−p̄:
   `n1 = [ zA·√(p̄q̄(1+1/r)) + zP·√(p1q1 + p2q2/r) ]² / (p1−p2)²`.
   **Unpooled**: both radicals `√(p1q1 + p2q2/r)`. **Arcsine (Cohen's h)** (r=1 only): `n/arm = 2(zA+zP)²/h²`.
   **Continuity correction** (Fleiss): r=1: `n' = (n/4)·(1+√(1+4/(n·|Δ|)))²`; general r: `n1' = (n1/4)·(1+√(1+2(r+1)/(r·n1·|Δ|)))²`.
   Power at given n: invert the same expressions (evaluate the two-term normal power: `power = Φ( (|Δ|·√n1 − zA·√(p̄q̄(1+1/r))) / √(p1q1+p2q2/r) )` for pooled; analogous unpooled/arcsine).
4. **Non-inferiority (means)** — one-sided by construction; margin m>0; true difference δ (positive = new arm better). `n1 = (z_{1−α}+zP)²·σ²·(1+1/r)/(δ+m)²`. Power: `Φ((δ+m)/SE − z_{1−α})`, SE=σ√((1+1/r)/n1). UI: direction picker + diagram (§10.7); tails forced 1; default α=0.025 with note.
5. **TOST equivalence (means)** — requires |δ| < m else return descriptive error. Power: `max(0, Φ((m−δ)/SE − z_{1−α}) + Φ((m+δ)/SE − z_{1−α}) − 1)`. Solve smallest n1 with power ≥ target.
6. **NI / TOST proportions** — same shifted-null structure with unpooled SE `√(p1q1 + p2q2/r)/√n1`.
7. **McNemar (paired proportions)** — Connor (1987): with Δ=|p2−p1| and discordant proportion `pd = p1q2 + p2q1 − 2ρ√(p1q1p2q2)`:
   `n_pairs = [ zA·√pd + zP·√(pd − Δ²) ]² / Δ²`.
   Validate ρ against `mcnemarRhoBounds` (all four cell probabilities p11=p1p2+Cov, p10=p1q2−Cov, p01=p2q1−Cov, p00=q1q2+Cov must be ≥0); out-of-range → descriptive error listing the admissible interval. Label: "large-sample approximation."
8. **k-arm pairwise** — m = comparisonsCount(k, scheme); α' = α/m; run calculator 1 or 3 at α'. The entered effect is the MINIMUM effect you care to detect on any tested pair (copy explains). Note: Bonferroni conservative; Dunnett would be slightly cheaper for vs-control (copy, not computed).
9. **Omnibus ANOVA** — Cohen's f; exact noncentral F. Derive helper: `anovaFFromMeans` (equal weights).
10. **Factorial a×b** — from cell-means grid + σ: grand mean μ̄; row effects αi (row marginal − μ̄); col effects βj; interaction γij = μij − μ̄ − αi − βj. `fA = √(Σαi²/a)/σ`, `fB = √(Σβj²/b)/σ`, `fAB = √(Σγij²/(ab))/σ`. Per-term solve via noncentral F (df: a−1, b−1, (a−1)(b−1); error df N−ab). Result = table of all three terms: derived f, min n per CELL, total N, achieved power. Highlight the interaction row.
11. **Omnibus χ² across k arms (binary DV)** — from per-arm projected proportions p_i (equal arm sizes): p̄ = mean(p_i); cells π1_i = p_i/k, π0_i = q_i/k under H1; under H0 p̄/k, q̄/k. `w = √(Σ_{2k cells} (π_H1 − π_H0)²/π_H0)`; df = k−1; λ = w²N; exact noncentral χ².
12. **Correlation** — Fisher z: `n = ceil( ((zA+zP)/atanh(|r|))² + 3 )`. Power at n inverts. Method label "Fisher z approximation"; note neighbor tools using exact-t agree at these magnitudes after ceiling.
13. **Precision mode** — 95% default confidence. Means: `n = (z·σ/h)²` then iterate with t: `n ← (tInv(1−(1−conf)/2, df(n))·σ/h)²` to fixed point (≤10 rounds), ceil. Mean difference per arm: same with σ√2. Proportion: `n = z²·p·q/h²` with "worst case p=0.5" toggle. Difference of proportions per arm: `n = z²·(p1q1+p2q2)/h²`. Copy: the width caveat — a CI planned this way has roughly a coin-flip chance of coming out wider; state it, don't hide it (§11.7).
14. **Derive helpers** — §4 list; each visible calculator embeds the relevant helper inline ("Don't know d? Type what you expect:").
15. **Adjustment pipeline (order is normative and displayed):** multiplicity adjusts α before solving → rank-inflation toggle (t-family only) → cluster DEFF → attrition → per-arm ceil. Result card shows the chain: base n → ×DEFF (m, ICC) → ÷(1−attrition) → final per-arm + clusters-per-arm when clustered.

## 7. Receipts (canonical self-test set)

Receipts run (a) in node via `test/receipts-runner.cjs` and (b) live in the page's Receipts panel via `runReceipts()`. Two classes:

- **PINNED** — values below are canon. If your engine disagrees, YOUR ENGINE IS WRONG (or the spec is — stop and report; never adjust a pinned value).
- **PIN-AT-VERIFICATION (PAV)** — compute with the finished engine, REPORT the value in your final report; the orchestrator independently verifies before it enters the shipped receipts table. Ship them marked `source: "engine, independently verified 2026-07-19"`.

Defaults unless stated: α=.05, two-tailed, power=.80, ratio=1.

| ID | Case | Expected | Source/method |
|----|------|----------|---------------|
| R01 | t2, d=0.5 → n/arm | **64** | pwr/G*Power canon (63.77) |
| R02 | t2, d=0.2 → n/arm | **394** | pwr (393.41) |
| R03 | t2, d=0.8 → n/arm | **26** | pwr (25.52) |
| R04 | t2, d=0.5, power=.90 → n/arm | **86** | pwr (85.03) |
| R05 | t2, d=0.5, one-tailed → n/arm | **51** | pwr (50.15) |
| R06 | paired, dz=0.5 → n pairs | **34** | pwr one-sample/paired (33.37) |
| R07 | dzFromPaired(0.5, 1, 0.5) = 0.5 → n = R06 | **34** | identity (σ_d=σ at ρ=.5) |
| R08 | one-sample, d=0.2 → n | **197** | pwr (196.22) |
| R09 | ANOVA k=3, f=0.25 → n/group | **53** | pwr (52.39) |
| R10 | ANOVA k=4, f=0.25 → n/group (N) | **45 (180)** | G*Power manual example |
| R11 | ANOVA k=2, f=0.25 ≡ t2 d=0.5 | **64/group** | identity F(1,ν)=t² |
| R12 | χ², w=0.3, df=2 → N | **108** | pwr (107.05) |
| R13 | χ², w=0.1, df=1 → N | **785** | pwr (784.90) |
| R14 | χ², w=0.5, df=1 → N | **32** | pwr (31.39) |
| R15 | correlation r=0.3 → n | **85** | Fisher z (84.93); exact-t neighbors agree at ceiling |
| R16 | props arcsine, .10 vs .20 → n/arm | **195** | h=0.28379; pwr.2p (194.91) |
| R17 | props arcsine, .50 vs .65 → n/arm | **170** | h=0.30468 (169.10) |
| R18 | props pooled-z, .10 vs .20 → n/arm | PAV (expect ≈199–200) | Fleiss formula |
| R19 | props pooled-z + CC, .10 vs .20 → n/arm | PAV (≈ R18 + ~19) | Fleiss CC |
| R20 | NI means, δ=0, m=0.4σ, α=.025 → n/arm | **99** | normal approx (98.11) |
| R21 | TOST means, δ=0, m=0.5σ, α=.05 → n/arm | **69** | normal approx (68.54) |
| R22 | McNemar p .50→.65, ρ=.3 → n pairs | PAV | Connor formula |
| R23 | cluster: base 64, m=20, ICC=.05 → n/arm; clusters/arm; warning | **125; 7; warning fires** | DEFF=1.95 |
| R24 | deff(m=1, any icc) = 1; applyCluster no-op | **exact** | identity |
| R25 | factorial 2×2 cells [[0,0],[0,0.5]], σ=1 → fA=fB=fAB=0.125; interaction n/cell | **f exact 0.125; n/cell PAV (expect ≈126)** | §6.10 derivation |
| R26 | pooled-z (no CC) vs χ² df=1 same n, .10 vs .20 | PAV equality | identity z²=χ² |
| R27 | mdeT2(64) round-trip | **d ∈ [0.495, 0.505]** | inversion consistency |
| R28 | precision mean-diff σ=1, h=0.2, 95% → n/arm | PAV (expect ≈194–197, exact-t) | §6.13 |
| R29 | powerT2(64, 0.5) | **∈ [0.8010, 0.8020]** | pwr reports 0.8015 |
| R30 | nctCdf sanity: nctCdf(1.6449, 1e6, 0) ≈ 0.95 | **±1e-4** | large-df → normal |

Rounding rule everywhere: per-arm/total n = `ceil` at the final step of each stage of the §6.15 pipeline (stages display their own ceils).

> **⚑ R08 CORRECTION — forward-only, pending Krystal ratification** (2026-07-19, Claudisegna A.-L. Frame Assayer v01 / Claude Opus 4.8, orchestrator, resuming the usage-limited build). The pinned R08 value **197** above is a **citation error** and the shipped receipts use **199**. Evidence, verified two independent ways: (1) the parenthetical "pwr (196.22)" is exactly the *asymptotic* one-sample n = (z₀.₉₇₅+z₀.₈₀)²/d² = 196.222 → ceil 197 — i.e. the author took the asymptotic ceiling, not the exact noncentral-t that an actual `pwr.t.test(type="one.sample")` computes; (2) the exact value is 198.14 → ceil **199**, reproduced *without the engine* by the Guenther small-sample correction 196.222 + z₀.₉₇₅²/2 = 198.143, and *with* the engine by the same `solveNT1` machinery that reproduces pinned R06 (d=0.5 → 34) exactly. The original **197** is preserved in the table above for provenance; the engine's `runReceipts()` R08 row carries the correction inline. **Revert = restore 197 in both places.** This is the only pinned value changed from the frozen spec.

## 8. Wizard (guided mode)

Four questions, plain language, big option cards. Routing is deterministic; every route states a one-sentence WHY ("Paired t: the same units appear in every condition, so each unit serves as its own control").

- **Q1 outcome:** "a number you'll average" / "a yes-or-no you'll count" / "ranks or ratings (clearly non-normal)" / "time-to-event, counts, or something else"
- **Q2 structure:** "2 conditions" / "3+ conditions, one factor" / "a grid of two factors" / "one group vs a known benchmark" / "two measurements on the same units (association)"
- **Q3 units:** "each unit experiences ONE condition" / "the same units experience EVERY condition" + independent checkbox "my units come in natural groups/clusters (classrooms, prompts, sites)"
- **Q4 goal:** "show a difference" / "show no-worse-than (needs a margin)" / "show equivalence (needs a margin)" / "estimate a quantity precisely"

Routing notes: ordinal/non-normal → route to the t-family with rank-inflation toggle ON + caveat, or escape hatch if severely non-normal; "something else" → Boundaries panel (escape hatch); equivalence/NI + factorial → advise reducing to the two-arm contrast of interest; cluster checkbox → cluster fields pre-armed on the destination calculator; association → correlation. Wizard hands off into the Calculator with fields pre-set and the route explanation pinned at top.

## 9. Planner (design canvas)

Sections (left-rail entries under "Plan"; all fields persist to the active plan):

1. **Question & hypotheses** — research question; hypotheses list `{text, direction: two-sided|greater|less, primary: bool}`. Exactly one primary enforced by UI nudge (multiplicity copy §11.5). Direction feeds tails.
2. **Variables** — IVs: `{name, levels[], between|within}`. DVs: `{name, type: continuous|binary|ordinal|other, instrument, primary}`. DV type feeds routing; ≥2 DVs triggers primary-endpoint nudge.
3. **Conditions grid** — auto-generated cross of IV levels (1 IV → arms list; 2 IVs → labeled a×b grid feeding §6.10; ≥3 IVs → grid shown but power routed to escape-hatch note "power the pre-registered contrasts or simulate"). Per cell: what differs / what's held constant.
4. **Assignment & randomization** — unit of randomization (text + "same as unit of analysis?" toggle → mismatch arms cluster fields + copy); method: simple / blocked / stratified-blocked; **working generator**: seed field (default 1234), N (prefill from calculator result), arm labels from grid, permuted block sizes (k and 2k mixed, toggle for fixed k), optional strata labels (one independent sequence per stratum); outputs preview table (first 20) + per-arm counts + CSV download (data URI). PRNG: mulberry32(seed) + Fisher-Yates. Deterministic: same inputs → same list.
5. **Confound register** — rows `{confound, mechanism, mitigation: hold-constant|randomize|measure-and-adjust|balance}`.
6. **Validity threats** — four checklists (construct / internal / external / statistical-conclusion), each item a prompt with a free-text response and a "addressed / accepted-risk / n/a" state. Construct items first (e.g., "What construct does the primary DV claim to measure? What would a skeptic say it actually measures?", "Does the manipulation move ONLY the construct you intend?"). Statistical-conclusion items link to the calculator ("Is the study powered for the SMALLEST effect you'd still care about?").
7. **Analysis plan** — auto-drafted from calculator state (test, α, tails, correction, adjustments) + free fields: exclusion rules, missing-data plan, secondary analyses.
8. **Power** — the calculator embedded, bound to the plan (k from grid, DV type from Variables, cluster fields from Assignment). Manual override allowed with a "plan ≠ calculator" flag.
9. **Export & save** — named plans (localStorage `eep.plans.v1`, autosave 500ms debounce, plan switcher: new/duplicate/rename/delete); JSON export/import (round-trip exact); **preregistration export**: one markdown document of the full plan + the power paragraph (§10.8 template) — copy button + `.md` download.

## 10. UI/UX

### 10.1 Layout
Left rail (collapsible): Wizard / Plan (9 sections) / Calculator / Receipts / Boundaries. Main canvas per selection. Header: title + unit-label field + plan switcher. Footer: version, receipts status chip (live count, links to panel), "made under the closed-form boundary" one-liner.

### 10.2 Design tokens
`--bg:#0a0d15; --panel:#121829; --panel2:#0f1422; --edge:#26304a; --ink:#e9eefb; --muted:#8b95ad; --cyan:#3ce6e6; --magenta:#ff5fc4; --lime:#a6ff5c; --amber:#ffd166; --violet:#b692ff; --orange:#ff9d5c;`
Glow: `box-shadow: 0 0 18px -4px <accent>`. Body text = ink; neon = accents/headings/data only. Section accents: Wizard=violet, Plan=cyan, Calculator=lime, Receipts=amber, Boundaries=magenta. Result numbers: big (clamp(2.2rem,4vw,3.2rem)) lime with glow.

### 10.3 Components
Option cards (wizard), field rows with inline explainer toggles ("?" buttons — §11 copy), slider+number pairs (α, power, effect, ρ, ICC, attrition), segmented chips (α: .05/.01/.10/custom; power: .80/.90/.95/custom; tails; method), result card (§10.6), collapsible "show the math" per calculator (renders the formula + the numbers substituted — plain text/Unicode, no LaTeX lib).

### 10.4 Power curve (SVG, no libs)
X = n/arm over [max(2, round(0.25·n*)), round(2.5·n*)], 60 sample points; Y = power 0–1. Curve in section accent; target-power horizontal dashed line; solution vertical line; hover crosshair with `(n, power)` readout; click sets n and flips card into "power at this n" mode. Second tab: MDE vs n (same interactions). Third tab when budget fields filled: affordable-N marker. Budget mini-fields: cost/unit, fixed cost, budget → affordable N → MDE at it.

### 10.5 Charts a11y
Curve readouts mirrored in a visually-hidden live region; all chart interactions keyboard-reachable (arrow keys move the crosshair).

### 10.6 Result card
Big per-arm n; total N; clusters/arm when clustered; achieved power at integer n; method line ("exact noncentral t" / "normal approximation (pooled z)" / "Connor large-sample approximation"); adjustment chain rendered as chips: `base 64 → ×1.95 DEFF → ÷0.9 attrition → 139/arm`; prereg sentence (§10.8) with copy button; assumption list link.

### 10.7 NI/TOST direction diagram
A number line showing the margin(s), null region, and where "better" lies, re-rendered from the direction picker. Static SVG, no animation.

### 10.8 Prereg power paragraph template
"To detect {effect phrase} on {DV name} with {power·100}% power at α = {α} ({tails}-tailed{, Bonferroni-adjusted to α′ = {α′} for {m} comparisons ({scheme})}), the {test name} requires {n} {unitLabel} per {arm|cell|group} ({N} total across {k} {arms|cells}). {Clustering: with {unitLabel} grouped in clusters of ~{m} (ICC = {icc}), the design effect is {deff}, requiring {nAdj} per arm = {clusters} clusters per arm.} {Attrition: anticipating {a}% attrition, we will recruit {nFinal} per arm ({NFinal} total).} Computed with exact noncentral-{t|F|χ²} methods ({method}); receipts panel documents validation."
Fill only applicable clauses; grammar must survive every combination (builder: test all branches).

## 11. Copy blocks (Builder B: adapt lightly, preserve content and stance)

1. **Boundary banner (Boundaries panel intro):** "Everything in this tool is closed-form and assumption-light. The moment your design outgrows these formulas, the honest tool stops pretending: simulate instead. Here's when, and how."
2. **SESOI doctrine (effect-size helpers):** "Power for the smallest effect you'd still care about, not the effect you hope for. Sources, best to worst: your pilot data → published effects in YOUR paradigm → Cohen's benchmarks (last resort — they were offered reluctantly and travel badly across fields)."
3. **4×/16× expandable (factorial):** "A 2×2 interaction is a difference of differences, so its estimate is noisier: at the same total N its standard error is twice a main effect's. Consequences — if the interaction contrast is as large as the marginal difference, you need ~4× the total N of the simple two-arm study that inspired it; if it's half as large (the common 'effect here, muted there' pattern), ~16×. Enter your predicted cell means above and the table computes the honest number for every term — no heuristic needed." Include the algebra in a nested details block (SE of difference-of-differences = 2σ·√(4/N)... show: Var(interaction contrast) = 4σ²·(4/N) vs Var(main contrast) = σ²·(4/N)).
4. **Cluster warning (< 8 clusters/arm):** "With this few clusters per arm, these formulas get optimistic: cluster-level noise dominates and the effective df are small. Consider more clusters (not more units per cluster — see the DEFF curve flatten), a cluster-level analysis, small-sample df corrections, or simulation."
5. **Multiplicity:** "Every comparison you'll report is another chance at a false positive. Bonferroni divides α across the m comparisons you commit to; it's conservative — the true cost is slightly lower (Dunnett for vs-control) — so treat these ns as safely sufficient. The cheapest fix is fewer confirmatory comparisons: one primary."
6. **Post-hoc power (Boundaries, refused):** "This tool won't compute 'observed power' from your results. Observed power is a deterministic transformation of your p-value — it adds no information, and 'we were underpowered' as an explanation of a null is circular (Hoenig & Heisey, 2001). Power is a DESIGN quantity: compute it before, against the smallest effect of interest."
7. **Precision caveat:** "A CI planned to half-width h has roughly a 50% chance of coming out somewhat wider (the width is itself random). If you need a guarantee, plan for a slightly smaller h."
8. **Escape hatch (Boundaries):** list of out-of-scope designs with WHY and WHERE: fully crossed random effects ("the same prompt set through every model, runs nested in both — power depends on variance components no form field should invent; simulate: R `simr`, `faux`, `Superpower`"), survival/time-to-event, counts, full ordinal proportional-odds, sequential/alpha-spending, Bayesian assurance. Close with the simulation principle: "simulate your analysis on data from the world you believe in; power = the fraction of simulated worlds where you detect it."
9. **Assumptions note (t-family):** equal-variance planning; Welch at analysis time costs almost nothing at near-equal arms; normality matters least at large n and most at tiny n (rank toggle exists; severe non-normality → simulate).
10. **One-tailed nudge:** "One-tailed is legitimate when direction is theory-committed AND pre-registered; it is not a discount code. If a reversed effect would still matter, stay two-tailed."
11. **McNemar label:** "Large-sample approximation from your projected marginals and correlation; the derived discordant probabilities are shown — sanity-check them."
12. **Rank toggle:** "+~5% n (ARE of rank tests vs t under shift alternatives ≈ 0.955). Varies with distribution shape; for heavy tails ranks can WIN — this toggle is a planning cushion, not a law."
13. **ICC guidance:** "Get ICC from pilots or literature in your setting. Education commonly 0.05–0.25 (classrooms); repeated prompts on a model can run far higher. When unsure, bracket: compute at low and high plausible ICC and look at the spread."

Tone: plain, direct, second person, zero hedging-as-decoration; label approximations plainly; Standard American English.

## 12. Acceptance criteria (both builders self-check before reporting)

1. `node test/receipts-runner.cjs` exits 0; all PINNED receipts pass; PAV values computed and listed in the report.
2. `node build.mjs` produces `index.html`; opened via `file://`: zero console errors, zero network requests, Receipts panel green.
3. Wizard: every path reaches a calculator or the Boundaries panel; each route shows its why-sentence.
4. Planner: 2-IV plan auto-builds the grid; grid k flows into calculator; randomizer deterministic (same seed → identical CSV twice); JSON export→import round-trips; prereg export grammatical in: base, +multiplicity, +cluster, +attrition, all-combined.
5. A11y/motion: keyboard-only pass of wizard→calculator→result; `prefers-reduced-motion` kills transitions; NO rotating/spinning/spiral anything; focus visible throughout; results announced via aria-live.
6. Every stat input has a working explainer; every result card shows method + adjustment chain + prereg sentence.
7. No library, no CDN, no web font, no external image (inline SVG only).

## 13. Division of labor and reporting

- **Builder A** — `src/stats-engine.js` + `test/receipts-runner.cjs`. Zero DOM. Runner: table output (id, desc, expected, got, pass), exit 1 on any pinned failure. Iterate until green. REPORT: PAV values, any spec discrepancies, function-by-function status.
- **Builder B** — `src/app-shell.html` complete per §§8–11 against the §4 API (engine absent until assembly — do NOT copy engine code or write a stub INTO the repo; a throwaway stub in the scratchpad for dev preview is fine). REPORT: acceptance-criteria status, any API needs beyond §4 (report, don't invent).
- **Orchestrator (Fable)** — build.mjs assembly, independent receipt verification, PAV ratification, ASAE gate, commits, delivery.
- Neither builder commits, pushes, or edits the other's files or docs/.

## Appendix A — worked constants

z(.975)=1.959964; z(.95)=1.644854; z(.90)=1.281552; z(.80)=0.841621.
(zA+zP)² at α=.05 2t, power .80: (1.959964+0.841621)² = 7.848886.
h(.10,.20)=0.283794; h(.50,.65)=0.304685.
atanh(0.3)=0.309520.
DEFF(m=20, ICC=.05)=1.95.
Factorial R25 grid [[0,0],[0,0.5]]: μ̄=0.125; row effects ±0.125; col effects ±0.125; γ = +0.125,−0.125,−0.125,+0.125 pattern; fA=fB=fAB=0.125 exactly.
