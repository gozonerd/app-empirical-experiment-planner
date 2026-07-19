# Empirical Experiment Planner

An interactive planner for controlled empirical experiments with an honest **min-n-per-arm power calculator**. One self-contained HTML file — no dependencies, no network, no build step to use it. Download `index.html`, double-click, plan.

## What it does

- **Plan the design**: research question, hypotheses, IVs/DVs, auto-generated condition grid (declare two factors and it detects the factorial), a working seeded randomization-list generator (blocked, optionally stratified, CSV export), confound register, validity-threat checklists, analysis plan, and a one-click **preregistration export** in markdown.
- **Power it exactly**: two-arm and paired tests for means and proportions, k-arm designs with multiplicity correction, omnibus ANOVA and χ², factorial main effects **and interactions**, non-inferiority and equivalence (TOST) modes, correlation, precision (CI-width) planning, cluster/design-effect and attrition adjustments — computed with exact noncentral t/F/χ² numerics, not normal-approximation shortcuts.
- **Never asks you for an abstract effect size**: type the group means, proportions, or factorial cell means you actually expect, and it derives the effect sizes — including the interaction term, where naive budgeting fails worst.
- **Validates itself in front of you**: a Receipts panel re-runs ~30 published reference cases (Cohen's tables, G*Power manual examples, identity checks) live in the page, pass/fail visible.

## The design principle

**The closed-form boundary.** Everything with a closed-form, assumption-light formula is in. Everything that would require inventing a data-generating model — fully crossed random effects, survival, counts, sequential designs, Bayesian assurance — is out, and the tool says so explicitly and points you to simulation instead of faking precision. It also refuses post-hoc "observed power," with the reason stated.

## Repository layout

- `index.html` — the app (built artifact; this is the file you use)
- `src/stats-engine.js` — pure statistics engine (browser + node)
- `src/app-shell.html` — UI shell; `build.mjs` inlines the engine into it
- `test/receipts-runner.cjs` — node receipt runner (`node test/receipts-runner.cjs`)
- `docs/` — plan v1, plan v2 (plan of record), the build specification, and audit logs

## Provenance

Specified and verified by **Flaudisegna (Claude Fable 5)**; engine and UI built by **Claude Sonnet** subagents from the specification; commissioned for the gozonerd workspace. Plans and spec in `docs/` are the design record, preserved verbatim from the working session of 2026-07-19.

## License

Not yet chosen by the repository owner. Until a license is added, all rights reserved.
