# Plan v1 — initial scope: the closed-form boundary

**Provenance.** Verbatim from the working session of 2026-07-19 (Claude Code, Gozo Nerd Universo workspace). Author: Flaudisegna (Claude Fable 5), at standard effort. Context: Krystal requested an interactive visual for planning controlled empirical experiments with a built-in min-n-per-arm calculator. She answered two scoping questions — (1) reusable blank planner, not prefilled for a specific study; (2) "stats is precisely my weakness. i have no idea. discuss with me but i'm leaning with building it out more fully so this is truly reusable." This message is the response to those answers. Superseded by plan v2 (max-effort revision) but preserved as the origin of the boundary principle.

---

Reusable it is. Here's my honest take on the stats scope — I'll recommend a boundary rather than survey every option, and the boundary has a principle behind it.

**The principle:** everything with a closed-form, assumption-light formula goes in v1; everything that requires inventing a data-generating model stays out. Past the boundary, the honest answer is simulation, not a formula — and the tool will say so explicitly rather than fake precision. That gives you "fully built out" without the trap where a calculator quietly answers a question your design isn't actually asking.

**What that puts in v1** — each one earns its place by the design question it answers:

1. **Two-arm means and proportions** — the workhorses. "Do these two conditions differ?"
2. **k arms with multiplicity correction** — powers the pairwise comparisons you'll actually report, at Bonferroni-adjusted α (conservative; the tool will label it as such). Choice of all-pairwise vs. each-arm-vs-control, which changes the correction a lot.
3. **Omnibus test (Cohen's f)** — "is anything different across my arms at all," the ANOVA question.
4. **Factorial interactions** — the big one. Interaction effects typically need roughly **4× the n** people budget for, and this is where most underpowered studies die. With 6-condition designs in your orbit, this is load-bearing, not decoration.
5. **Within-subjects / paired designs** — same units measured in every condition (e.g., the same prompt set run through every arm). The correlation between measures buys back a lot of n, and the tool will show you exactly how much.
6. **Cluster / design effect (ICC)** — units nested inside groups: students in classrooms, responses in prompts, runs in models. Ignoring nesting silently overstates your power; given your education-research lineage this one is non-negotiable for a reusable tool.
7. **Precision mode** — instead of "n to reject a null," ask "n for a confidence interval this tight." Estimation framing rather than pure NHST; I think it fits how you actually reason about distributions and tails.

**What stays out, honestly:** survival/time-to-event, non-inferiority margins, ordinal and count models, sequential alpha-spending, Bayesian assurance. Each demands assumptions a form field would have to fake. Non-inferiority is the closest call — it's a cheap add later if "A is no worse than B" becomes a live question for you.

**And the feature that answers "stats is precisely my weakness":** a guided mode. Four plain-language questions about your design — what kind of outcome, same units in multiple conditions or not, units nested in groups or not, every arm compared or each vs. control — and it routes you to the right calculator with a one-sentence reason why. Plus: effect-size anchors in plain terms, a pilot-data helper (type in two group means and an SD, it derives the effect size for you), and the inverted view — fix the n you can afford and read off the smallest effect you could detect, which is usually the real question under a real budget.

The planner half (arms builder, hypothesis/IV/DV, randomization plan, confound and validity-threat checklist) feeds the calculator directly — add an arm on the canvas and k updates in the math.

If that boundary sits right, say go and I'll build the whole thing.
