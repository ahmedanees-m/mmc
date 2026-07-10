# Abstract

MMC is a mechanistic discovery loop for gene-regulatory circuits. It reads the Zhu 2025
genome-scale CD4+ T-cell Perturb-seq atlas, proposes interpretable, runnable models of
regulatory circuits (signed edges and simulatable logic) through a reasoning step, fits them to
the measured perturbation responses, and evaluates each against simple baselines on held-out
data. The reasoning step sets structure and logic; the optimizer sets magnitudes.

On the Th2 / GATA3 axis the loop reconstructs a coherent circuit (in-sample Pearson 0.93). On
the cytokine-production module it proposes STK11 / LKB1 as a chemokine regulator; an
edge-ablation control confirms the STK11 edges are individually supported, comparable to
textbook edges, but the model built from them does not beat a linear baseline on held-out data
(DE-overlap 0.18 versus 0.45, separated CIs), so the hypothesis is reported as proposed but not
certified. Across a powered corpus of 25 proposals and 76 distinct novel hypotheses over nine
runs and two conditions, none is in a model that beats a linear baseline held-out (Wilson 95% CI
[0, 4.8%]). A reasoning-versus-search comparison shows the reasoning step is not equivalent to
random structure search, producing better-fitting structure than a random structure of equal
edge count on the informative modules (mean in-sample Pearson 0.20 versus 0.07). Edge-level
support does not imply predictive advantage over a simple baseline; the module-level held-out
evaluation is what distinguishes them.

MMC includes a limit map: the regime boundary of where mechanistic and AI models beat simple
baselines on held-out data. On single-knockdown steady-state data no held-out advantage is
observed in any measured regime; the remaining source of advantage is non-additivity, and a
combinatorial test (Norman K562) shows that a model fit on single perturbations does not recover
it either.

Scope: CD4+ T cells, these modules, this atlas. No prediction win and no disease discovery is
claimed. The contribution is a method for producing interpretable, testable regulatory
hypotheses and a reproducible boundary for when mechanistic and AI models are reliable.
