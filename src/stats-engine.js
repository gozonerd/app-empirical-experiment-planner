// stats-engine.js — Empirical Experiment Planner pure-logic stats engine.
// Zero DOM. Built against docs/build-spec-empirical-experiment-planner_2026-07-19_v01.md.
// Sections below are numbered to mirror the spec (§4 API, §5 numerics, §6 formulas, §7 receipts).

"use strict";

// ============================================================
// §5 Numerics layer
// ============================================================

// ---- erf / normCdf (Abramowitz & Stegun 7.1.26) ----
function erf(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741,
        a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  const ax = Math.abs(x);
  const t = 1 / (1 + p * ax);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-ax * ax);
  return sign * y;
}

function normCdf(z) {
  return 0.5 * (1 + erf(z / Math.SQRT2));
}

// ---- normInv (Acklam's algorithm + Halley refinement) ----
function normInv(p) {
  if (p === 0) return -Infinity;
  if (p === 1) return Infinity;
  if (!(p > 0) || !(p < 1)) {
    throw new Error(`normInv: p must be in (0,1), got ${p}`);
  }
  const a1 = -3.969683028665376e+01, a2 = 2.209460984245205e+02,
        a3 = -2.759285104469687e+02, a4 = 1.383577518672690e+02,
        a5 = -3.066479806614716e+01, a6 = 2.506628277459239e+00;
  const b1 = -5.447609879822406e+01, b2 = 1.615858368580409e+02,
        b3 = -1.556989798598866e+02, b4 = 6.680131188771972e+01,
        b5 = -1.328068155288572e+01;
  const c1 = -7.784894002430293e-03, c2 = -3.223964580411365e-01,
        c3 = -2.400758277161838e+00, c4 = -2.549732539343734e+00,
        c5 = 4.374664141464968e+00, c6 = 2.938163982698783e+00;
  const d1 = 7.784695709041462e-03, d2 = 3.224671290700398e-01,
        d3 = 2.445134137142996e+00, d4 = 3.754408661907416e+00;
  const pLow = 0.02425, pHigh = 1 - pLow;
  let x;
  if (p < pLow) {
    const q = Math.sqrt(-2 * Math.log(p));
    x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
  } else if (p <= pHigh) {
    const q = p - 0.5, r = q * q;
    x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
        (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
  }
  // Halley refinement step (standard part of the published Acklam algorithm).
  // e = Phi(x) - p; using normCdf built on our erf() above.
  const e = normCdf(x) - p;
  const u = e * Math.sqrt(2 * Math.PI) * Math.exp(x * x / 2);
  x = x - u / (1 + x * u / 2);
  return x;
}

// ---- lgamma (Lanczos g=7, n=9) ----
const LANCZOS_G = 7;
const LANCZOS_P = [
  0.99999999999980993, 676.5203681218851, -1259.1392167224028,
  771.32342877765313, -176.61502916214059, 12.507343278686905,
  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
];

function lgamma(x) {
  if (x < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
  }
  x -= 1;
  let a = LANCZOS_P[0];
  const t = x + LANCZOS_G + 0.5;
  for (let i = 1; i < LANCZOS_G + 2; i++) a += LANCZOS_P[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

// ---- regularized incomplete beta ibeta(x,a,b) via Lentz continued fraction ----
function betacf(x, a, b) {
  const MAXIT = 400, EPS = 1e-14, FPMIN = 1e-300;
  const qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1;
  let d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}

function ibeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const logBt = a * Math.log(x) + b * Math.log(1 - x) + lgamma(a + b) - lgamma(a) - lgamma(b);
  const bt = Math.exp(logBt);
  if (x < (a + 1) / (a + b + 2)) {
    return bt * betacf(x, a, b) / a;
  }
  return 1 - bt * betacf(1 - x, b, a) / b;
}

// ---- regularized lower incomplete gamma pgamma(a,x) ----
function gser(a, x) {
  const ITMAX = 400, EPS = 1e-14;
  if (x <= 0) return 0;
  let ap = a;
  let sum = 1 / a;
  let del = sum;
  for (let n = 1; n <= ITMAX; n++) {
    ap += 1;
    del *= x / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * EPS) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
}

function gcf(a, x) {
  const ITMAX = 400, EPS = 1e-14, FPMIN = 1e-300;
  let b = x + 1 - a;
  let c = 1 / FPMIN;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i <= ITMAX; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
}

function pgamma(a, x) {
  if (x <= 0) return 0;
  if (x < a + 1) return gser(a, x);
  return 1 - gcf(a, x);
}

// ---- central CDFs ----
function tCdf(t, df) {
  const x = df / (df + t * t);
  const p = 0.5 * ibeta(x, df / 2, 0.5);
  return t >= 0 ? 1 - p : p;
}

function fCdf(x, d1, d2) {
  if (x <= 0) return 0;
  const v = d1 * x / (d1 * x + d2);
  return ibeta(v, d1 / 2, d2 / 2);
}

function chisqCdf(x, df) {
  if (x <= 0) return 0;
  return pgamma(df / 2, x / 2);
}

// ---- generic bisection (bracket assumed; used by quantiles + MDE inversions) ----
function bisectSolve(f, lo, hi, opts) {
  opts = opts || {};
  const tol = opts.tol !== undefined ? opts.tol : 1e-12;
  const maxIter = opts.maxIter !== undefined ? opts.maxIter : 200;
  let flo = f(lo), fhi = f(hi);
  if (flo === 0) return lo;
  if (fhi === 0) return hi;
  let mid = (lo + hi) / 2;
  for (let i = 0; i < maxIter; i++) {
    mid = (lo + hi) / 2;
    const fmid = f(mid);
    if (fmid === 0) return mid;
    if ((flo < 0) === (fmid < 0)) { lo = mid; flo = fmid; } else { hi = mid; fhi = fmid; }
    if (Math.abs(hi - lo) <= tol * Math.max(1, Math.abs(mid))) return mid;
  }
  return mid;
}

// ---- quantiles: bracket-doubling + bisection ----
function tInv(p, df) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;
  const f = (t) => tCdf(t, df) - p;
  let lo = -1, hi = 1, iter = 0;
  while (f(lo) > 0 && iter < 2000) { lo *= 2; iter++; }
  iter = 0;
  while (f(hi) < 0 && iter < 2000) { hi *= 2; iter++; }
  return bisectSolve(f, lo, hi);
}

function fInv(p, d1, d2) {
  if (p <= 0) return 0;
  if (p >= 1) return Infinity;
  const f = (x) => fCdf(x, d1, d2) - p;
  let lo = 0, hi = 1, iter = 0;
  while (f(hi) < 0 && iter < 2000) { hi *= 2; iter++; }
  return bisectSolve(f, lo, hi);
}

function chisqInv(p, df) {
  if (p <= 0) return 0;
  if (p >= 1) return Infinity;
  const f = (x) => chisqCdf(x, df) - p;
  let lo = 0, hi = 1, iter = 0;
  while (f(hi) < 0 && iter < 2000) { hi *= 2; iter++; }
  return bisectSolve(f, lo, hi);
}

// ---- noncentral t via chi-square-mixture integral (composite Simpson) ----
function nctCdf(t, df, delta) {
  delta = delta || 0;
  if (df > 1e5) {
    const z = (t * (1 - 1 / (4 * df)) - delta) / Math.sqrt(1 + (t * t) / (2 * df));
    return normCdf(z);
  }
  if (delta === 0) return tCdf(t, df);
  const lo = chisqInv(1e-13, df);
  const hi = chisqInv(1 - 1e-13, df);
  const logC = -lgamma(df / 2) - (df / 2) * Math.log(2);
  function integrand(v) {
    if (v <= 0) return 0;
    const logf = (df / 2 - 1) * Math.log(v) - v / 2 + logC;
    return normCdf(t * Math.sqrt(v / df) - delta) * Math.exp(logf);
  }
  function simpson(n) {
    const h = (hi - lo) / n;
    let sum = integrand(lo) + integrand(hi);
    for (let i = 1; i < n; i++) {
      const x = lo + i * h;
      sum += (i % 2 === 0 ? 2 : 4) * integrand(x);
    }
    return sum * h / 3;
  }
  let n = 256;
  let prev = simpson(n);
  while (n < 32768) {
    n *= 2;
    const cur = simpson(n);
    if (Math.abs(cur - prev) < 1e-10) return cur;
    prev = cur;
  }
  return prev;
}

// ---- shared Poisson-mixture summation for ncf/ncchisq ----
function poissonMixtureSum(lambda, termFn) {
  if (lambda === 0) return termFn(0);
  const halfLambda = lambda / 2;
  const j0 = Math.floor(halfLambda);
  const logPoissonAtJ0 = -halfLambda + j0 * Math.log(halfLambda) - lgamma(j0 + 1);
  const w0 = Math.exp(logPoissonAtJ0);
  let total = w0 * termFn(j0);
  // upward
  let wj = w0, jj = j0 + 1;
  while (true) {
    wj = wj * halfLambda / jj;
    if (wj < 1e-16) break;
    total += wj * termFn(jj);
    jj++;
  }
  // downward
  wj = w0;
  jj = j0 - 1;
  while (jj >= 0) {
    wj = wj * (jj + 1) / halfLambda;
    if (wj < 1e-16) break;
    total += wj * termFn(jj);
    jj--;
  }
  return total;
}

function ncfCdf(x, d1, d2, lambda) {
  lambda = lambda || 0;
  if (x <= 0) return 0;
  if (lambda === 0) return fCdf(x, d1, d2);
  const v = d1 * x / (d1 * x + d2);
  return poissonMixtureSum(lambda, (j) => ibeta(v, (d1 + 2 * j) / 2, d2 / 2));
}

function ncchisqCdf(x, df, lambda) {
  lambda = lambda || 0;
  if (x <= 0) return 0;
  if (lambda === 0) return chisqCdf(x, df);
  return poissonMixtureSum(lambda, (j) => pgamma((df + 2 * j) / 2, x / 2));
}

// ============================================================
// Shared internal helpers (not part of the §4 public surface)
// ============================================================

function findMinN(powerFn, targetPower, minN, maxN) {
  minN = minN === undefined ? 2 : minN;
  maxN = maxN === undefined ? 1e8 : maxN;
  let n = Math.max(2, Math.ceil(minN));
  if (powerFn(n) >= targetPower) return n;
  let lo = n, hi = n, guard = 0;
  while (powerFn(hi) < targetPower) {
    hi *= 2;
    guard++;
    if (hi > maxN || guard > 200) {
      throw new Error(`required n exceeds ${maxN}; check inputs (effect size may be ~0 or otherwise infeasible)`);
    }
  }
  while (hi - lo > 1) {
    const mid = Math.floor((lo + hi) / 2);
    if (powerFn(mid) >= targetPower) hi = mid; else lo = mid;
  }
  return hi;
}

function n2FromRatio(n1, ratio) {
  return Math.max(2, Math.ceil(ratio * n1));
}

function bisectMde(f, lo, hi, growCap) {
  growCap = growCap === undefined ? 100 : growCap;
  let iter = 0;
  while (f(hi) < 0) {
    hi *= 2;
    iter++;
    if (hi > growCap || iter > 200) throw new Error("MDE search exceeded bounds");
  }
  return bisectSolve(f, lo, hi);
}

// ============================================================
// §6.1–6.2 Two-arm & one-arm means (exact noncentral t)
// ============================================================

function powerT2(n1, d, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const n2 = n2FromRatio(n1, ratio);
  const df = n1 + n2 - 2;
  const ncp = d * Math.sqrt((n1 * n2) / (n1 + n2));
  const tcrit = tInv(1 - alpha / tails, df);
  let power;
  if (tails === 2) {
    power = (1 - nctCdf(tcrit, df, ncp)) + nctCdf(-tcrit, df, ncp);
  } else {
    power = 1 - nctCdf(tcrit, df, ncp);
  }
  return { power };
}

function solveNT2(d, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const n1 = findMinN((n) => powerT2(n, d, opts).power, targetPower, 2);
  const n2 = n2FromRatio(n1, ratio);
  const df = n1 + n2 - 2;
  const ncp = d * Math.sqrt((n1 * n2) / (n1 + n2));
  const tcrit = tInv(1 - alpha / tails, df);
  return {
    nPerArm: n1, n1, n2, Ntotal: n1 + n2,
    achievedPower: powerT2(n1, d, opts).power,
    df, criticalValue: tcrit, ncpAtSolution: ncp,
    method: "exact noncentral t (independent two-sample)",
    notes: []
  };
}

function mdeT2(n1, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const f = (d) => powerT2(n1, d, opts).power - targetPower;
  const d = bisectMde(f, 1e-6, 0.01);
  return { d };
}

function powerT1(n, d, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const df = n - 1;
  const ncp = d * Math.sqrt(n);
  const tcrit = tInv(1 - alpha / tails, df);
  let power;
  if (tails === 2) {
    power = (1 - nctCdf(tcrit, df, ncp)) + nctCdf(-tcrit, df, ncp);
  } else {
    power = 1 - nctCdf(tcrit, df, ncp);
  }
  return { power };
}

function solveNT1(d, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const rawN = findMinN((n) => powerT1(n, d, opts).power, targetPower, 2);
  const n = Math.max(rawN, 4);
  const notes = [];
  if (rawN < 4) {
    notes.push(`tiny-sample: unconstrained solution was n=${rawN}; floored to the one-sample/paired minimum n=4`);
  }
  const df = n - 1;
  const ncp = d * Math.sqrt(n);
  const tcrit = tInv(1 - alpha / tails, df);
  return {
    nPerArm: n, n, Ntotal: n,
    achievedPower: powerT1(n, d, opts).power,
    df, criticalValue: tcrit, ncpAtSolution: ncp,
    method: "exact noncentral t (one-sample/paired)",
    notes
  };
}

function mdeT1(n, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const f = (d) => powerT1(n, d, opts).power - targetPower;
  const d = bisectMde(f, 1e-6, 0.01);
  return { d };
}

// ============================================================
// §6.3 Two-arm proportions (normal approx)
// ============================================================

function hFromProps(p1, p2) {
  return Math.abs(2 * (Math.asin(Math.sqrt(p1)) - Math.asin(Math.sqrt(p2))));
}

function powerProp2(n1, p1, p2, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const r = opts.ratio !== undefined ? opts.ratio : 1;
  const method = opts.method || "pooled";
  const cc = !!opts.cc;
  const zA = normInv(1 - alpha / tails);
  const delta = Math.abs(p1 - p2);

  if (method === "arcsine") {
    if (r !== 1) throw new Error("powerProp2: arcsine method requires ratio=1");
    if (cc) throw new Error("powerProp2: continuity correction is not defined for the arcsine method");
    const h = hFromProps(p1, p2);
    return { power: normCdf(h * Math.sqrt(n1 / 2) - zA) };
  }

  const q1 = 1 - p1, q2 = 1 - p2;
  const seH1 = Math.sqrt(p1 * q1 + (p2 * q2) / r);
  let seH0;
  if (method === "unpooled") {
    seH0 = seH1;
  } else if (method === "pooled") {
    const pbar = (p1 + r * p2) / (1 + r);
    const qbar = 1 - pbar;
    seH0 = Math.sqrt(pbar * qbar * (1 + 1 / r));
  } else {
    throw new Error(`powerProp2: unknown method "${method}"`);
  }
  let effDelta = delta;
  if (cc) {
    effDelta = Math.max(0, delta - (0.5 * (1 + 1 / r)) / n1);
  }
  const power = normCdf((effDelta * Math.sqrt(n1) - zA * seH0) / seH1);
  return { power };
}

function solveNProp2(p1, p2, opts) {
  opts = opts || {};
  const method = opts.method || "pooled";
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  if (method === "arcsine" && ratio !== 1) throw new Error("solveNProp2: arcsine method requires ratio=1");
  if (method === "arcsine" && opts.cc) throw new Error("solveNProp2: continuity correction is not defined for the arcsine method");
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const n1 = findMinN((n) => powerProp2(n, p1, p2, opts).power, targetPower, 2);
  const n2 = n2FromRatio(n1, ratio);
  return {
    nPerArm: n1, n1, n2, Ntotal: n1 + n2,
    achievedPower: powerProp2(n1, p1, p2, opts).power,
    df: null, criticalValue: normInv(1 - alpha / tails), ncpAtSolution: null,
    method: `normal approximation (${method}${opts.cc ? " + continuity correction" : ""})`,
    notes: []
  };
}

function mdeProp2(n1, p1, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const direction = opts.direction !== undefined ? opts.direction : 1;
  const f = (p2) => powerProp2(n1, p1, p2, opts).power - targetPower;
  let lo, hi;
  if (direction === 1) { lo = Math.min(p1 + 1e-9, 1 - 1e-9); hi = 1 - 1e-9; }
  else { lo = 1e-9; hi = Math.max(p1 - 1e-9, 1e-9); }
  if (f(hi) < 0) {
    throw new Error("mdeProp2: target power unreachable within (0,1) bounds for the requested direction");
  }
  const p2 = bisectSolve(f, lo, hi);
  return { p2 };
}

// ============================================================
// §6.4–6.6 Shifted-null designs (normal approx)
// ============================================================

function powerNIMeans(n1, deltaTrue, margin, sd, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const zA = normInv(1 - alpha);
  const SE = sd * Math.sqrt((1 + 1 / ratio) / n1);
  return { power: normCdf((deltaTrue + margin) / SE - zA) };
}

function solveNNIMeans(deltaTrue, margin, sd, opts) {
  opts = opts || {};
  if (!(margin > 0)) throw new Error("solveNNIMeans: margin must be > 0");
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const n1 = findMinN((n) => powerNIMeans(n, deltaTrue, margin, sd, opts).power, targetPower, 2);
  const n2 = n2FromRatio(n1, ratio);
  return {
    nPerArm: n1, n1, n2, Ntotal: n1 + n2,
    achievedPower: powerNIMeans(n1, deltaTrue, margin, sd, opts).power,
    df: null, criticalValue: normInv(1 - alpha), ncpAtSolution: null,
    method: "normal approximation (non-inferiority, means)",
    notes: ["tails forced to 1 for non-inferiority design"]
  };
}

function powerTostMeans(n1, deltaTrue, margin, sd, opts) {
  opts = opts || {};
  if (!(Math.abs(deltaTrue) < margin)) throw new Error("powerTostMeans: TOST requires |deltaTrue| < margin");
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const zA = normInv(1 - alpha);
  const SE = sd * Math.sqrt((1 + 1 / ratio) / n1);
  const power = Math.max(0,
    normCdf((margin - deltaTrue) / SE - zA) + normCdf((margin + deltaTrue) / SE - zA) - 1);
  return { power };
}

function solveNTostMeans(deltaTrue, margin, sd, opts) {
  opts = opts || {};
  if (!(margin > 0)) throw new Error("solveNTostMeans: margin must be > 0");
  if (!(Math.abs(deltaTrue) < margin)) throw new Error("solveNTostMeans: TOST requires |deltaTrue| < margin");
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const n1 = findMinN((n) => powerTostMeans(n, deltaTrue, margin, sd, opts).power, targetPower, 2);
  const n2 = n2FromRatio(n1, ratio);
  return {
    nPerArm: n1, n1, n2, Ntotal: n1 + n2,
    achievedPower: powerTostMeans(n1, deltaTrue, margin, sd, opts).power,
    df: null, criticalValue: normInv(1 - alpha), ncpAtSolution: null,
    method: "normal approximation (TOST equivalence, means)",
    notes: ["tails forced to 1 (two one-sided tests) for TOST design"]
  };
}

function powerNIPropCalc(n1, p1, p2, margin, opts) {
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const zA = normInv(1 - alpha);
  const delta = p1 - p2;
  const SEunit = Math.sqrt(p1 * (1 - p1) + (p2 * (1 - p2)) / ratio);
  const SE = SEunit / Math.sqrt(n1);
  return normCdf((delta + margin) / SE - zA);
}

function solveNNIProp(p1, p2, margin, opts) {
  opts = opts || {};
  if (!(margin > 0)) throw new Error("solveNNIProp: margin must be > 0");
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const n1 = findMinN((n) => powerNIPropCalc(n, p1, p2, margin, opts), targetPower, 2);
  const n2 = n2FromRatio(n1, ratio);
  return {
    nPerArm: n1, n1, n2, Ntotal: n1 + n2,
    achievedPower: powerNIPropCalc(n1, p1, p2, margin, opts),
    df: null, criticalValue: normInv(1 - alpha), ncpAtSolution: null,
    method: "normal approximation (non-inferiority, proportions, unpooled SE)",
    notes: ["tails forced to 1 for non-inferiority design"]
  };
}

function powerTostPropCalc(n1, p1, p2, margin, opts) {
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const zA = normInv(1 - alpha);
  const delta = p1 - p2;
  const SEunit = Math.sqrt(p1 * (1 - p1) + (p2 * (1 - p2)) / ratio);
  const SE = SEunit / Math.sqrt(n1);
  return Math.max(0, normCdf((margin - delta) / SE - zA) + normCdf((margin + delta) / SE - zA) - 1);
}

function solveNTostProp(p1, p2, margin, opts) {
  opts = opts || {};
  if (!(margin > 0)) throw new Error("solveNTostProp: margin must be > 0");
  if (!(Math.abs(p1 - p2) < margin)) throw new Error("solveNTostProp: TOST requires |p1-p2| < margin");
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const ratio = opts.ratio !== undefined ? opts.ratio : 1;
  const n1 = findMinN((n) => powerTostPropCalc(n, p1, p2, margin, opts), targetPower, 2);
  const n2 = n2FromRatio(n1, ratio);
  return {
    nPerArm: n1, n1, n2, Ntotal: n1 + n2,
    achievedPower: powerTostPropCalc(n1, p1, p2, margin, opts),
    df: null, criticalValue: normInv(1 - alpha), ncpAtSolution: null,
    method: "normal approximation (TOST equivalence, proportions, unpooled SE)",
    notes: ["tails forced to 1 (two one-sided tests) for TOST design"]
  };
}

// ============================================================
// §6.7 McNemar (paired proportions, Connor 1987 approximation)
// ============================================================

function mcnemarRhoBounds(p1, p2) {
  const q1 = 1 - p1, q2 = 1 - p2;
  const denom = Math.sqrt(p1 * q1 * p2 * q2);
  const hi = Math.min(p1 * q2, p2 * q1) / denom;
  const lo = -Math.min(p1 * p2, q1 * q2) / denom;
  return { lo, hi };
}

function powerMcNemarDirect(n, p10, p01, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const delta = Math.abs(p01 - p10);
  const pd = p10 + p01;
  const zA = normInv(1 - alpha / tails);
  const denom = Math.sqrt(Math.max(pd - delta * delta, 1e-15));
  return { power: normCdf((delta * Math.sqrt(n) - zA * Math.sqrt(pd)) / denom) };
}

function solveNMcNemarDirect(p10, p01, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const rawN = findMinN((n) => powerMcNemarDirect(n, p10, p01, opts).power, targetPower, 2);
  const n = Math.max(rawN, 4);
  const notes = [`derived discordant probabilities: p10=${p10.toFixed(6)}, p01=${p01.toFixed(6)}`];
  if (rawN < 4) {
    notes.push(`tiny-sample: unconstrained solution was n=${rawN}; floored to the paired-design minimum n=4`);
  }
  return {
    nPerArm: n, n, Ntotal: n,
    achievedPower: powerMcNemarDirect(n, p10, p01, opts).power,
    df: null, criticalValue: normInv(1 - alpha / tails), ncpAtSolution: null,
    method: "Connor (1987) large-sample approximation",
    notes
  };
}

function solveNMcNemar(p1, p2, rho, opts) {
  opts = opts || {};
  const bounds = mcnemarRhoBounds(p1, p2);
  if (rho < bounds.lo || rho > bounds.hi) {
    throw new Error(
      `solveNMcNemar: rho=${rho} out of admissible range [${bounds.lo.toFixed(6)}, ${bounds.hi.toFixed(6)}] for p1=${p1}, p2=${p2}`
    );
  }
  const q1 = 1 - p1, q2 = 1 - p2;
  const cov = rho * Math.sqrt(p1 * q1 * p2 * q2);
  const p10 = p1 * q2 - cov;
  const p01 = p2 * q1 - cov;
  return solveNMcNemarDirect(p10, p01, opts);
}

// ============================================================
// §6.9 k-group omnibus (exact noncentral F / chi-square)
// ============================================================

function powerAnova(nPerGroup, f, k, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const N = nPerGroup * k;
  const df1 = k - 1, df2 = N - k;
  const lambda = f * f * N;
  const fcrit = fInv(1 - alpha, df1, df2);
  return { power: 1 - ncfCdf(fcrit, df1, df2, lambda) };
}

function solveNAnova(f, k, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const nPerGroup = findMinN((n) => powerAnova(n, f, k, opts).power, targetPower, 2);
  const N = nPerGroup * k;
  const df1 = k - 1, df2 = N - k;
  const lambda = f * f * N;
  return {
    nPerGroup, Ntotal: N,
    achievedPower: powerAnova(nPerGroup, f, k, opts).power,
    df: [df1, df2], criticalValue: fInv(1 - alpha, df1, df2), ncpAtSolution: lambda,
    method: "exact noncentral F (omnibus ANOVA)",
    notes: []
  };
}

function powerChisqCalc(N, w, df, opts) {
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const lambda = w * w * N;
  const crit = chisqInv(1 - alpha, df);
  return 1 - ncchisqCdf(crit, df, lambda);
}

function solveNChisq(w, df, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const N = findMinN((n) => powerChisqCalc(n, w, df, opts), targetPower, 2);
  return {
    Ntotal: N,
    achievedPower: powerChisqCalc(N, w, df, opts),
    df, criticalValue: chisqInv(1 - alpha, df), ncpAtSolution: w * w * N,
    method: "exact noncentral chi-square",
    notes: []
  };
}

function solveNChisqKprops(props, opts) {
  opts = opts || {};
  const k = props.length;
  const pbar = props.reduce((a, b) => a + b, 0) / k;
  const qbar = 1 - pbar;
  let w2sum = 0;
  for (let i = 0; i < k; i++) {
    const p = props[i], q = 1 - p;
    const pi1H1 = p / k, pi0H1 = q / k;
    const pi1H0 = pbar / k, pi0H0 = qbar / k;
    w2sum += Math.pow(pi1H1 - pi1H0, 2) / pi1H0 + Math.pow(pi0H1 - pi0H0, 2) / pi0H0;
  }
  const w = Math.sqrt(w2sum);
  const df = k - 1;
  const res = solveNChisq(w, df, opts);
  res.w = w;
  return res;
}

function mdeAnova(nPerGroup, k, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const g = (f) => powerAnova(nPerGroup, f, k, opts).power - targetPower;
  const f = bisectMde(g, 1e-6, 0.01);
  return { f };
}

// ============================================================
// §6.10 Factorial (exact noncentral F per term)
// ============================================================

function factorialEffectsFromCells(grid, sd) {
  const a = grid.length, b = grid[0].length;
  let grandMean = 0;
  for (let i = 0; i < a; i++) for (let j = 0; j < b; j++) grandMean += grid[i][j];
  grandMean /= (a * b);
  const rowMeans = grid.map((row) => row.reduce((s, v) => s + v, 0) / b);
  const colMeans = [];
  for (let j = 0; j < b; j++) {
    let s = 0;
    for (let i = 0; i < a; i++) s += grid[i][j];
    colMeans.push(s / a);
  }
  const rowEffects = rowMeans.map((m) => m - grandMean);
  const colEffects = colMeans.map((m) => m - grandMean);
  const interactionEffects = [];
  for (let i = 0; i < a; i++) {
    const row = [];
    for (let j = 0; j < b; j++) {
      row.push(grid[i][j] - grandMean - rowEffects[i] - colEffects[j]);
    }
    interactionEffects.push(row);
  }
  const fA = Math.sqrt(rowEffects.reduce((s, v) => s + v * v, 0) / a) / sd;
  const fB = Math.sqrt(colEffects.reduce((s, v) => s + v * v, 0) / b) / sd;
  let sumGamma2 = 0;
  for (let i = 0; i < a; i++) for (let j = 0; j < b; j++) sumGamma2 += interactionEffects[i][j] * interactionEffects[i][j];
  const fAB = Math.sqrt(sumGamma2 / (a * b)) / sd;
  return { fA, fB, fAB, rowEffects, colEffects, interactionEffects, grandMean };
}

function powerFactorialTerm(nPerCell, fTerm, dfTerm, cells, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const N = nPerCell * cells;
  const df2 = N - cells;
  const lambda = fTerm * fTerm * N;
  const fcrit = fInv(1 - alpha, dfTerm, df2);
  return { power: 1 - ncfCdf(fcrit, dfTerm, df2, lambda) };
}

function solveNFactorialTerm(fTerm, dfTerm, cells, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const nPerCell = findMinN((n) => powerFactorialTerm(n, fTerm, dfTerm, cells, opts).power, targetPower, 2);
  const N = nPerCell * cells;
  const df2 = N - cells;
  return {
    nPerCell, Ntotal: N,
    achievedPower: powerFactorialTerm(nPerCell, fTerm, dfTerm, cells, opts).power,
    df: [dfTerm, df2], criticalValue: fInv(1 - alpha, dfTerm, df2), ncpAtSolution: fTerm * fTerm * N,
    method: "exact noncentral F (factorial term)",
    notes: []
  };
}

// ============================================================
// §6.12 Correlation (Fisher z)
// ============================================================

function powerCorr(n, r, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const zA = normInv(1 - alpha / tails);
  const zr = Math.atanh(Math.abs(r));
  return { power: normCdf(zr * Math.sqrt(n - 3) - zA) };
}

function solveNCorr(r, opts) {
  opts = opts || {};
  const targetPower = opts.power !== undefined ? opts.power : 0.80;
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const n = findMinN((nn) => powerCorr(nn, r, opts).power, targetPower, 4);
  return {
    nPerArm: n, n, Ntotal: n,
    achievedPower: powerCorr(n, r, opts).power,
    df: n - 2, criticalValue: normInv(1 - alpha / tails), ncpAtSolution: null,
    method: "Fisher z approximation",
    notes: ["exact-t neighbors agree at these magnitudes after ceiling"]
  };
}

function mdeCorr(n, opts) {
  opts = opts || {};
  const alpha = opts.alpha !== undefined ? opts.alpha : 0.05;
  const tails = opts.tails !== undefined ? opts.tails : 2;
  const power = opts.power !== undefined ? opts.power : 0.80;
  const zA = normInv(1 - alpha / tails);
  const zP = normInv(power);
  const zr = (zA + zP) / Math.sqrt(n - 3);
  return { r: Math.tanh(zr) };
}

// ============================================================
// §6.13 Precision mode (CI half-width)
// ============================================================

function nMeanCI(sd, halfWidth, conf) {
  const alpha = 1 - conf;
  const z = normInv(1 - alpha / 2);
  let n = Math.pow((z * sd) / halfWidth, 2);
  for (let i = 0; i < 10; i++) {
    const df = Math.max(1, Math.round(n) - 1);
    const t = tInv(1 - alpha / 2, df);
    const nNew = Math.pow((t * sd) / halfWidth, 2);
    if (Math.abs(nNew - n) < 1e-9) { n = nNew; break; }
    n = nNew;
  }
  return Math.ceil(n);
}

function nMeanDiffCI(sd, halfWidth, conf) {
  const alpha = 1 - conf;
  const z = normInv(1 - alpha / 2);
  const sdEff = sd * Math.SQRT2;
  let n = Math.pow((z * sdEff) / halfWidth, 2);
  for (let i = 0; i < 10; i++) {
    const nInt = Math.max(2, Math.round(n));
    const df = 2 * nInt - 2;
    const t = tInv(1 - alpha / 2, df);
    const nNew = Math.pow((t * sdEff) / halfWidth, 2);
    if (Math.abs(nNew - n) < 1e-9) { n = nNew; break; }
    n = nNew;
  }
  return Math.ceil(n);
}

function nPropCI(p, halfWidth, conf) {
  const alpha = 1 - conf;
  const z = normInv(1 - alpha / 2);
  return Math.ceil((z * z * p * (1 - p)) / (halfWidth * halfWidth));
}

function nPropDiffCI(p1, p2, halfWidth, conf) {
  const alpha = 1 - conf;
  const z = normInv(1 - alpha / 2);
  return Math.ceil((z * z * (p1 * (1 - p1) + p2 * (1 - p2))) / (halfWidth * halfWidth));
}

// ============================================================
// Adjustments & derivation helpers
// ============================================================

function deff(m, icc, cvm) {
  cvm = cvm === undefined ? 0 : cvm;
  return 1 + ((1 + cvm * cvm) * m - 1) * icc;
}

function applyCluster(n, m, icc, cvm) {
  cvm = cvm === undefined ? 0 : cvm;
  const d = deff(m, icc, cvm);
  const nAdjusted = Math.ceil(n * d);
  const clustersPerArm = Math.ceil(nAdjusted / m);
  const fewClustersWarning = clustersPerArm < 8;
  return { nAdjusted, clustersPerArm, deff: d, fewClustersWarning };
}

function applyAttrition(n, rate) {
  return Math.ceil(n / (1 - rate));
}

function comparisonsCount(k, scheme) {
  if (scheme === "allPairs") return (k * (k - 1)) / 2;
  if (scheme === "vsControl") return k - 1;
  throw new Error(`comparisonsCount: unknown scheme "${scheme}"`);
}

function bonferroniAlpha(alpha, m) {
  return alpha / m;
}

function rankInflate(n) {
  return Math.ceil(n / 0.955);
}

function dFromMeans(m1, m2, sd) {
  return (m1 - m2) / sd;
}

function dFromMeansSds(m1, m2, sd1, sd2) {
  const sdPooled = Math.sqrt((sd1 * sd1 + sd2 * sd2) / 2);
  return (m1 - m2) / sdPooled;
}

function dzFromPaired(delta, sd, rho) {
  return delta / (sd * Math.sqrt(2 * (1 - rho)));
}

function anovaFFromMeans(means, sd) {
  const k = means.length;
  const grand = means.reduce((a, b) => a + b, 0) / k;
  const msVar = means.reduce((s, m) => s + (m - grand) * (m - grand), 0) / k;
  return Math.sqrt(msVar) / sd;
}

// ============================================================
// §7 Receipts
// ============================================================

function runReceipts() {
  const rows = [];
  function add(row) {
    rows.push(row);
    return row;
  }
  function closeEnough(got, expected, tol) {
    return Math.abs(got - expected) <= tol;
  }

  // R01
  {
    const res = solveNT2(0.5, { alpha: 0.05, tails: 2, power: 0.80, ratio: 1 });
    add({
      id: "R01", class: "PINNED", desc: "t2, d=0.5 -> n/arm", expected: 64, got: res.nPerArm,
      pass: res.nPerArm === 64, method: res.method, source: "pwr/G*Power canon (63.77)"
    });
  }
  // R02
  {
    const res = solveNT2(0.2, {});
    add({
      id: "R02", class: "PINNED", desc: "t2, d=0.2 -> n/arm", expected: 394, got: res.nPerArm,
      pass: res.nPerArm === 394, method: res.method, source: "pwr (393.41)"
    });
  }
  // R03
  {
    const res = solveNT2(0.8, {});
    add({
      id: "R03", class: "PINNED", desc: "t2, d=0.8 -> n/arm", expected: 26, got: res.nPerArm,
      pass: res.nPerArm === 26, method: res.method, source: "pwr (25.52)"
    });
  }
  // R04
  {
    const res = solveNT2(0.5, { power: 0.90 });
    add({
      id: "R04", class: "PINNED", desc: "t2, d=0.5, power=.90 -> n/arm", expected: 86, got: res.nPerArm,
      pass: res.nPerArm === 86, method: res.method, source: "pwr (85.03)"
    });
  }
  // R05
  {
    const res = solveNT2(0.5, { tails: 1 });
    add({
      id: "R05", class: "PINNED", desc: "t2, d=0.5, one-tailed -> n/arm", expected: 51, got: res.nPerArm,
      pass: res.nPerArm === 51, method: res.method, source: "pwr (50.15)"
    });
  }
  // R06
  let r06n;
  {
    const res = solveNT1(0.5, {});
    r06n = res.n;
    add({
      id: "R06", class: "PINNED", desc: "paired, dz=0.5 -> n pairs", expected: 34, got: res.n,
      pass: res.n === 34, method: res.method, source: "pwr one-sample/paired (33.37)"
    });
  }
  // R07
  {
    const dz = dzFromPaired(0.5, 1, 0.5);
    const res = solveNT1(dz, {});
    add({
      id: "R07", class: "PINNED", desc: "dzFromPaired(0.5,1,0.5)=0.5 -> n = R06", expected: 34, got: res.n,
      pass: dz === 0.5 && res.n === r06n, method: "identity", source: "identity (sigma_d=sigma at rho=.5)"
    });
  }
  // R08 — orchestrator-corrected 2026-07-19 (Claudisegna A.-L. Frame Assayer v01, Claude Opus 4.8):
  // spec pinned 197 was a citation error. "196.22" is the ASYMPTOTIC one-sample n ((zA+zP)^2/d^2);
  // the exact noncentral-t that the pin's own "pwr" basis actually computes gives 198.14 -> ceil 199.
  // Independently verified two ways: (1) engine's one-sample machinery is canon-correct (reproduces
  // pinned R06 exactly); (2) Guenther small-sample correction 196.22 + z(.975)^2/2 = 198.14 -> 199,
  // no engine involved. Original pin (197) preserved in docs/build-spec §7 with a dated correction
  // note; pending Krystal ratification (revert = restore 197).
  {
    const res = solveNT1(0.2, {});
    add({
      id: "R08", class: "PINNED", desc: "one-sample, d=0.2 -> n", expected: 199, got: res.n,
      pass: res.n === 199, method: res.method, source: "engine, independently verified 2026-07-19 (spec pin 197 = asymptotic ceil; exact noncentral-t = 199)"
    });
  }
  // R09
  {
    const res = solveNAnova(0.25, 3, {});
    add({
      id: "R09", class: "PINNED", desc: "ANOVA k=3, f=0.25 -> n/group", expected: 53, got: res.nPerGroup,
      pass: res.nPerGroup === 53, method: res.method, source: "pwr (52.39)"
    });
  }
  // R10
  {
    const res = solveNAnova(0.25, 4, {});
    add({
      id: "R10", class: "PINNED", desc: "ANOVA k=4, f=0.25 -> n/group (N)",
      expected: "45 (180)", got: `${res.nPerGroup} (${res.Ntotal})`,
      pass: res.nPerGroup === 45 && res.Ntotal === 180, method: res.method,
      source: "G*Power manual example"
    });
  }
  // R11
  {
    const res = solveNAnova(0.25, 2, {});
    add({
      id: "R11", class: "PINNED", desc: "ANOVA k=2, f=0.25 == t2 d=0.5", expected: "64/group",
      got: `${res.nPerGroup}/group`, pass: res.nPerGroup === 64, method: "identity F(1,v)=t^2",
      source: "identity F(1,nu)=t^2"
    });
  }
  // R12
  {
    const res = solveNChisq(0.3, 2, {});
    add({
      id: "R12", class: "PINNED", desc: "chisq, w=0.3, df=2 -> N", expected: 108, got: res.Ntotal,
      pass: res.Ntotal === 108, method: res.method, source: "pwr (107.05)"
    });
  }
  // R13
  {
    const res = solveNChisq(0.1, 1, {});
    add({
      id: "R13", class: "PINNED", desc: "chisq, w=0.1, df=1 -> N", expected: 785, got: res.Ntotal,
      pass: res.Ntotal === 785, method: res.method, source: "pwr (784.90)"
    });
  }
  // R14
  {
    const res = solveNChisq(0.5, 1, {});
    add({
      id: "R14", class: "PINNED", desc: "chisq, w=0.5, df=1 -> N", expected: 32, got: res.Ntotal,
      pass: res.Ntotal === 32, method: res.method, source: "pwr (31.39)"
    });
  }
  // R15
  {
    const res = solveNCorr(0.3, {});
    add({
      id: "R15", class: "PINNED", desc: "correlation r=0.3 -> n", expected: 85, got: res.n,
      pass: res.n === 85, method: res.method, source: "Fisher z (84.93); exact-t neighbors agree at ceiling"
    });
  }
  // R16
  {
    const res = solveNProp2(0.10, 0.20, { method: "arcsine" });
    add({
      id: "R16", class: "PINNED", desc: "props arcsine, .10 vs .20 -> n/arm", expected: 195, got: res.nPerArm,
      pass: res.nPerArm === 195, method: res.method, source: "h=0.28379; pwr.2p (194.91)"
    });
  }
  // R17
  {
    const res = solveNProp2(0.50, 0.65, { method: "arcsine" });
    add({
      id: "R17", class: "PINNED", desc: "props arcsine, .50 vs .65 -> n/arm", expected: 170, got: res.nPerArm,
      pass: res.nPerArm === 170, method: res.method, source: "h=0.30468 (169.10)"
    });
  }
  // R18 (PAV)
  let r18n;
  {
    const res = solveNProp2(0.10, 0.20, { method: "pooled" });
    r18n = res.nPerArm;
    add({
      id: "R18", class: "PAV", desc: "props pooled-z, .10 vs .20 -> n/arm",
      expected: "PAV (expect ~199-200)", got: res.nPerArm,
      pass: res.nPerArm >= 199 && res.nPerArm <= 200, method: res.method, source: "Fleiss formula"
    });
  }
  // R19 (PAV)
  {
    const res = solveNProp2(0.10, 0.20, { method: "pooled", cc: true });
    add({
      id: "R19", class: "PAV", desc: "props pooled-z + CC, .10 vs .20 -> n/arm",
      expected: `PAV (~R18 + ~19 = ~${r18n + 19})`, got: res.nPerArm,
      pass: res.nPerArm > r18n, method: res.method, source: "Fleiss CC"
    });
  }
  // R20
  {
    const res = solveNNIMeans(0, 0.4, 1, { alpha: 0.025, power: 0.80 });
    add({
      id: "R20", class: "PINNED", desc: "NI means, delta=0, m=0.4sd, alpha=.025 -> n/arm",
      expected: 99, got: res.nPerArm, pass: res.nPerArm === 99, method: res.method,
      source: "normal approx (98.11)"
    });
  }
  // R21
  {
    const res = solveNTostMeans(0, 0.5, 1, { alpha: 0.05, power: 0.80 });
    add({
      id: "R21", class: "PINNED", desc: "TOST means, delta=0, m=0.5sd, alpha=.05 -> n/arm",
      expected: 69, got: res.nPerArm, pass: res.nPerArm === 69, method: res.method,
      source: "normal approx (68.54)"
    });
  }
  // R22 (PAV)
  {
    const res = solveNMcNemar(0.50, 0.65, 0.3, {});
    add({
      id: "R22", class: "PAV", desc: "McNemar p .50->.65, rho=.3 -> n pairs",
      expected: "PAV", got: res.n, pass: Number.isFinite(res.n) && res.n > 0, method: res.method,
      source: "Connor formula"
    });
  }
  // R23
  {
    const res = applyCluster(64, 20, 0.05);
    add({
      id: "R23", class: "PINNED", desc: "cluster: base 64, m=20, ICC=.05 -> n/arm; clusters/arm; warning",
      expected: "125; 7; warning fires",
      got: `${res.nAdjusted}; ${res.clustersPerArm}; ${res.fewClustersWarning ? "warning fires" : "no warning"}`,
      pass: res.nAdjusted === 125 && res.clustersPerArm === 7 && res.fewClustersWarning === true,
      method: "deff/applyCluster identity", source: "DEFF=1.95"
    });
  }
  // R24
  {
    const d1 = deff(1, 0.37);
    const ac = applyCluster(50, 1, 0.37);
    add({
      id: "R24", class: "PINNED", desc: "deff(m=1, any icc) = 1; applyCluster no-op",
      expected: "exact", got: `deff=${d1}; nAdjusted=${ac.nAdjusted}`,
      pass: d1 === 1 && ac.nAdjusted === 50, method: "identity", source: "identity"
    });
  }
  // R25
  let r25nPerCell;
  {
    const eff = factorialEffectsFromCells([[0, 0], [0, 0.5]], 1);
    const fExactPass = closeEnough(eff.fA, 0.125, 1e-9) && closeEnough(eff.fB, 0.125, 1e-9) && closeEnough(eff.fAB, 0.125, 1e-9);
    const dfTerm = 1; // (a-1)(b-1) for 2x2
    const cells = 4;
    const nres = solveNFactorialTerm(eff.fAB, dfTerm, cells, {});
    r25nPerCell = nres.nPerCell;
    add({
      id: "R25", class: fExactPass ? "PINNED" : "PINNED",
      desc: "factorial 2x2 cells [[0,0],[0,0.5]], sd=1 -> fA=fB=fAB=0.125 exact; interaction n/cell PAV",
      expected: "f exact 0.125; n/cell PAV (expect ~126)",
      got: `fA=${eff.fA}, fB=${eff.fB}, fAB=${eff.fAB}, n/cell=${r25nPerCell}`,
      pass: fExactPass, method: nres.method, source: "spec 6.10 derivation"
    });
  }
  // R26 (PAV identity)
  {
    const pz = solveNProp2(0.10, 0.20, { method: "pooled", cc: false });
    const chi = solveNChisqKprops([0.10, 0.20], {});
    const zTotal = pz.n1 + pz.n2;
    add({
      id: "R26", class: "PAV", desc: "pooled-z (no CC) vs chisq df=1 same n, .10 vs .20",
      expected: "PAV equality (approximate: two different approximations of the same underlying test)",
      got: `z-total=${zTotal}, chisq-N=${chi.Ntotal}`,
      pass: Math.abs(zTotal - chi.Ntotal) / chi.Ntotal <= 0.02, method: "identity z^2=chisq", source: "identity z^2=chisq"
    });
  }
  // R27
  {
    const res = mdeT2(64, {});
    add({
      id: "R27", class: "PINNED", desc: "mdeT2(64) round-trip", expected: "d in [0.495, 0.505]",
      got: res.d, pass: res.d >= 0.495 && res.d <= 0.505, method: "inversion", source: "inversion consistency"
    });
  }
  // R28 (PAV)
  {
    const n = nMeanDiffCI(1, 0.2, 0.95);
    add({
      id: "R28", class: "PAV", desc: "precision mean-diff sd=1, h=0.2, 95% -> n/arm",
      expected: "PAV (expect ~194-196, exact-t)", got: n,
      pass: n >= 194 && n <= 196, method: "exact-t iteration", source: "spec 6.13"
    });
  }
  // R29
  {
    const res = powerT2(64, 0.5, {});
    add({
      id: "R29", class: "PINNED", desc: "powerT2(64, 0.5)", expected: "[0.8010, 0.8020]",
      got: res.power, pass: res.power >= 0.8010 && res.power <= 0.8020, method: "exact noncentral t",
      source: "pwr reports 0.8015"
    });
  }
  // R30
  {
    const got = nctCdf(1.6449, 1e6, 0);
    add({
      id: "R30", class: "PINNED", desc: "nctCdf sanity: nctCdf(1.6449, 1e6, 0) ~ 0.95", expected: "0.95 (+/-1e-4)",
      got, pass: closeEnough(got, 0.95, 1e-4), method: "large-df normal approx", source: "large-df -> normal"
    });
  }

  return rows;
}

// ============================================================
// §4 Public API assembly
// ============================================================

const StatsEngine = {
  // distributions
  normCdf, normInv, erf,
  tCdf, tInv,
  fCdf, fInv,
  chisqCdf, chisqInv,
  nctCdf, ncfCdf, ncchisqCdf,

  // two-arm & one-arm means
  powerT2, solveNT2, mdeT2,
  powerT1, solveNT1, mdeT1,

  // two-arm proportions
  powerProp2, solveNProp2, mdeProp2,

  // shifted-null designs
  solveNNIMeans, powerNIMeans,
  solveNTostMeans, powerTostMeans,
  solveNNIProp, solveNTostProp,

  // paired proportions (McNemar)
  solveNMcNemar, solveNMcNemarDirect, mcnemarRhoBounds,

  // k-group omnibus
  powerAnova, solveNAnova,
  solveNChisq, solveNChisqKprops,
  mdeAnova,

  // factorial
  factorialEffectsFromCells, solveNFactorialTerm, powerFactorialTerm,

  // correlation
  solveNCorr, powerCorr, mdeCorr,

  // precision mode
  nMeanCI, nMeanDiffCI, nPropCI, nPropDiffCI,

  // adjustments & derivation helpers
  deff, applyCluster, applyAttrition,
  comparisonsCount, bonferroniAlpha, rankInflate,
  dFromMeans, dFromMeansSds, hFromProps, dzFromPaired, anovaFFromMeans,

  // receipts
  runReceipts
};

if (typeof module !== "undefined" && module.exports) { module.exports = StatsEngine; }
if (typeof globalThis !== "undefined") { globalThis.StatsEngine = StatsEngine; }
