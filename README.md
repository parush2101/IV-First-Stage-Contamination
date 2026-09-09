# What the IV First Stage Cannot Detect

Replication code and manuscripts for the project *What the IV First Stage Cannot
Detect*, by Parush Arora (Department of Economics, Ashoka University).

In just-identified linear IV with a binary instrument, a binary treatment, and
covariates entered linearly, an exact decomposition writes the 2SLS bias as a
covariance between the curvature of the instrument propensity and the covariate
level functions. The same object builds the first-stage covariance the
Sanderson–Windmeijer conditional *F* is computed from, so a single nuisance both
biases the estimand and inflates the reported strength. No first-stage diagnostic
is informative about the bias, which is pinned down only by the conditional law of
the outcome. The repository provides a directed conditional-moment test, a
bias-corrected estimator, a reportable contamination share, and the applications.

## Repository layout

```
paper/     LaTeX sources and compiled PDFs of both versions
code/      R scripts that generate every table and figure
figures/   figure outputs embedded in the long version
```

## Manuscripts (`paper/`)

| File | Version | Target |
|------|---------|--------|
| `contamination_bias_IV_AERinsights.{tex,pdf}` | Short | *AER: Insights* |
| `contamination_bias_IV.{tex,pdf}`             | Long  | *Journal of Econometrics* |

Both are self-contained: the bibliography is embedded via `thebibliography`, so no
`.bib` file is needed. Compile with `pdflatex` **twice**. The short version draws
its figures natively with `pgfplots`; the long version embeds the four PNGs in
`figures/` (referenced as `../figures/`), so compile it from inside `paper/`.

## Code (`code/`)

| Script | Produces |
|--------|----------|
| `iv_collinearity_mc.R`      | Collinearity Monte Carlo (Table 1, collinearity figure) |
| `iv_contamination_test.R`   | The directed conditional-moment test |
| `iv_contamination_share.R`  | Contamination share ϱ (figure) |
| `iv_headtohead.R`           | RESET vs. directed tests across designs (head-to-head) |
| `iv_svx_collinearity.R`     | Sanderson–Windmeijer conditional *F* under collinearity |
| `iv_svx_harmless.R`         | Harmless-curvature design |
| `iv_power_curve.R`          | Power curve (figure) |
| `iv_401k.R`                 | Headline 401(k) design |
| `iv_401k_robustCI.R`        | 401(k) Anderson–Rubin / weak-IV-robust interval |
| `iv_empirical_final.R`      | Applied panel (Card, Mroz, Cigarettes, AJR, Maimonides) |
| `iv_maimonides.R`           | Angrist–Lavy class-size design |

## Data

Datasets are **not** included in this repository. The scripts read their inputs by
filename and are documented in the code. Public sources:

- Wooldridge teaching datasets (`401ksubs`, `card`, `catholic`, `crime4`,
  `fertil2`, `openness`, `wage2`, `htv`) — available with Wooldridge, *Econometric
  Analysis of Cross Section and Panel Data*.
- `AER` R package built-ins (`PSID1976` / Mroz, `CigarettesSW`, `AJR`) — installed
  with the package.
- Angrist–Lavy Maimonides class-size data — MIT/Angrist data archive.

Point R at the directory holding the `.dta` files (`setwd()`), or edit the paths at
the top of each script.

## Author

Parush Arora, Ashoka University · parush.arora@ashoka.edu.in
