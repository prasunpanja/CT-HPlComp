---
title: "User Manual for CT-HPlComp"
author: "Prasun Panja"
---

# 1. Overview

This document is a simple user manual for the R code provided in `hpl(1).R`.

The code implements functions for testing pleiotropy between a genotype and two traits. The main procedures included are:

- `hoprr()` — permutation-based HoP-RR test.
- `plafiu()` — forward-regression based Pl-AFIU-style calculation.
- `marg_plei()` — marginal pleiotropy test with Bonferroni correction.
- Supporting numerical routines for finite-difference derivatives and Newton-Raphson estimation.
- A toy simulation used to generate genotype, marker, and trait data.

The main workflow is:

```text
Genotype + two traits
        |
        v
   hoprr() / plafiu()
        |
        v
 Test statistic / p-value / decision
```

# 2. Requirements

The supplied code uses base R functions and does not explicitly require an external R package.

A reasonably recent version of R is recommended.

To check the R version:

```{r}
R.version.string
```

# 3. Input Data Structure

The main functions expect:

- `gt`: genotype vector, generally coded as `0`, `1`, or `2`.
- `y`: a matrix/data frame containing two quantitative traits.
- `cv`: optional covariates.

For two traits, the expected structure of `y` is:

```text
        Trait 1   Trait 2
Sample 1    ...       ...
Sample 2    ...       ...
Sample 3    ...       ...
...
```

The genotype vector must contain the same number of observations as the rows of `y`.

For example:

```{r, eval=FALSE}
gt <- c(0, 1, 2, 1, 0)
y <- cbind(
  Trait1 = c(2.1, 3.4, 4.2, 2.8, 1.9),
  Trait2 = c(5.1, 6.2, 7.0, 5.7, 4.8)
)
```

# 4. Main Functions

## 4.1 `hoprr()`

The `hoprr()` function performs the HoP-RR procedure implemented in the supplied code.

Its basic syntax is:

```{r, eval=FALSE}
hoprr(gt, y, cv = NULL, n_perm = 1000, seed = 123)
```

### Arguments

| Argument | Description |
|---|---|
| `gt` | Genotype vector, coded as 0/1/2. |
| `y` | Matrix/data frame containing the two traits. |
| `cv` | Optional covariates. Default is `NULL`. |
| `n_perm` | Number of valid permutations. Default is 1000. |
| `seed` | Random-number seed. Default is 123. |

### What the function does

The function first fits the observed binomial reverse-regression model using `binomial_lrt()`.

It then:

1. Obtains the observed likelihood-ratio statistic.
2. Estimates genotype probabilities under the selected null model.
3. Generates genotype vectors under the null model using `rbinom()`.
4. Repeats the likelihood-ratio calculation for the permuted data.
5. Retains valid permutation results.
6. Calculates the empirical p-value.
7. Classifies the result as `"Significant"` when `p <= 0.05`.

The returned object contains:

```text
sample_size
est
lrt
p_value
decision
valid_permutations
total_attempted
```

### Example

```{r, eval=FALSE}
result <- hoprr(
  gt = Y[, 2],
  y = cbind(Y[, 3], Y[, 4]),
  n_perm = 1000,
  seed = 123
)

result
```

To extract individual results:

```{r, eval=FALSE}
result$p_value
result$lrt
result$decision
```

# 5. `binomial_lrt()`

`binomial_lrt()` is the core likelihood-ratio calculation used by `hoprr()`.

Syntax:

```{r, eval=FALSE}
binomial_lrt(gt, y, cv = NULL)
```

The function compares models in which the genotype is explained by:

- Trait 1 plus covariates,
- Trait 2 plus covariates,

and a full model containing both traits plus covariates.

The stronger null model is selected according to its likelihood, and the likelihood-ratio statistic is calculated against the full model.

The function returns:

```text
lrt
null_prob
null1
null2
full
```

This function is mainly a supporting function and normally does not need to be called directly by a user.

# 6. `lm_binom()`

`lm_binom()` fits the binomial regression used by the likelihood calculations.

Syntax:

```{r, eval=FALSE}
lm_binom(y, x, trials = 2)
```

The function:

1. Adds an intercept to the design matrix.
2. Defines the binomial-logistic score function.
3. Estimates the regression coefficients using Newton-Raphson.
4. Calculates fitted genotype probabilities.
5. Calculates the binomial log-likelihood.

It returns:

```text
estimate
prob
likelihood
```

This is an internal computational function.

# 7. Numerical Optimization Functions

## 7.1 `numeric_gradient()`

This function calculates a numerical gradient using finite differences.

```{r, eval=FALSE}
numeric_gradient(f, x, h1 = 1e-8)
```

It is a supporting numerical routine.

## 7.2 `jacobian()`

This function calculates the Jacobian matrix of a vector-valued function using finite differences.

```{r, eval=FALSE}
jacobian(f, x, h1 = 1e-8)
```

## 7.3 `slfn()`

`slfn()` solves a system of nonlinear equations using the Newton-Raphson method.

```{r, eval=FALSE}
slfn(f, x0, tol = 1e-6, max_iter = 100)
```

Arguments:

- `f`: vector-valued function.
- `x0`: initial parameter vector.
- `tol`: convergence tolerance.
- `max_iter`: maximum number of iterations.

If the Jacobian is singular or ill-conditioned, the function stops with an error.

If convergence is not reached within `max_iter`, the function also stops.

# 8. `plafiu()`

The `plafiu()` function performs two forward linear regressions:

```text
Trait 1 ~ Genotype + Trait 2 + covariates

Trait 2 ~ Genotype + Trait 1 + covariates
```

Syntax:

```{r, eval=FALSE}
plafiu(gt, y, cv = NULL)
```

The function returns:

```text
effect1
effect2
sd1
p1
p2
maxp
```

The `maxp` value is the maximum of the two regression p-values.

Example:

```{r, eval=FALSE}
plafiu(
  gt = Y[, 2],
  y = cbind(Y[, 3], Y[, 4])
)
```

# 9. `marg_plei()`

`marg_plei()` performs marginal genotype-trait tests and applies a Bonferroni threshold of 0.025 to each of the two tests.

Syntax:

```{r, eval=FALSE}
marg_plei(gt, y, cv)
```

The function compares:

```text
Genotype ~ covariates

versus

Genotype ~ Trait 1 + covariates
```

and

```text
Genotype ~ covariates

versus

Genotype ~ Trait 2 + covariates
```

The returned object contains:

```text
marginal1
marginal2
p1
p2
est1
est2
margD
```

The result is labelled `"Significant"` only when both marginal p-values are below 0.025.

# 10. Toy Simulation

The supplied script contains a complete toy-data simulation.

The main simulation parameters include:

```{r, eval=FALSE}
p1 <- 0.2
p2 <- 0.2
v1 <- 25
v2 <- 25
sm_size <- 500
```

Here:

- `p1` is the MAF of the trait locus.
- `p2` is the MAF of the marker locus.
- `v1` and `v2` specify the total variances of the two traits.
- `sm_size` is the sample size.

The simulation generates:

- a trait-locus genotype `A`,
- two quantitative traits,
- dichotomized versions of the traits,
- a marker genotype `M`.

The final simulated data matrix is:

```{r, eval=FALSE}
Y <- cbind(A, M, y[,1], y[,2])
```

Thus, the columns of `Y` are:

| Column | Meaning |
|---|---|
| `Y[,1]` | Trait locus genotype |
| `Y[,2]` | Marker genotype |
| `Y[,3]` | Quantitative Trait 1 |
| `Y[,4]` | Quantitative Trait 2 |

# 11. Running the Supplied Example

After running the complete script, the supplied implementation calls:

```{r, eval=FALSE}
hoprr(Y[,2], cbind(Y[,3], Y[,4]))
```

and

```{r, eval=FALSE}
plafiu(Y[,2], cbind(Y[,3], Y[,4]))
```

These calls use the simulated marker genotype as `gt` and the two simulated quantitative traits as `y`.

# 12. Applying the Functions to Your Own Data

Suppose your data frame is called `dat` and contains:

```text
SNP    Trait1    Trait2    Age    Sex
```

You can run:

```{r, eval=FALSE}
gt <- dat$SNP

y <- cbind(
  Trait1 = dat$Trait1,
  Trait2 = dat$Trait2
)

cv <- cbind(
  Age = dat$Age,
  Sex = dat$Sex
)
```

Then run HoP-RR:

```{r, eval=FALSE}
res <- hoprr(
  gt = gt,
  y = y,
  cv = cv,
  n_perm = 1000,
  seed = 123
)

res
```

For a larger permutation analysis, increase `n_perm`, for example:

```{r, eval=FALSE}
res <- hoprr(
  gt = gt,
  y = y,
  cv = cv,
  n_perm = 10000,
  seed = 123
)
```

# 13. Interpreting the HoP-RR Output

The most important fields are:

### `lrt`

The observed likelihood-ratio test statistic.

### `p_value`

The empirical permutation p-value calculated as the proportion of valid permuted statistics that exceed the observed statistic.

### `decision`

The supplied code uses:

```text
p_value <= 0.05  -> Significant
p_value > 0.05   -> Insignificant
```

### `valid_permutations`

The number of successful permutations retained. This should equal the requested `n_perm`.

### `total_attempted`

The total number of permutation attempts required to obtain the requested number of valid permutations.

# 14. Recommended Basic Workflow

For a single SNP:

```{r, eval=FALSE}
# 1. Prepare genotype
gt <- dat$SNP

# 2. Prepare the two traits
y <- cbind(dat$Trait1, dat$Trait2)

# 3. Prepare covariates if required
cv <- cbind(dat$Age, dat$Sex)

# 4. Run HoP-RR
res <- hoprr(
  gt = gt,
  y = y,
  cv = cv,
  n_perm = 1000,
  seed = 123
)

# 5. Inspect result
print(res)

# 6. Extract p-value
res$p_value
```

# 15. Important Practical Checks

Before running the analysis, check that:

1. `gt`, `y`, and `cv` have compatible sample sizes.
2. Genotypes are coded consistently as `0`, `1`, and `2`.
3. Trait values are numeric.
4. Covariates are appropriately encoded.
5. Missing values are handled before calling the functions.
6. The number of permutations is sufficiently large for the desired p-value resolution.
7. The random seed is kept fixed when reproducibility is required.

For example:

```{r, eval=FALSE}
length(gt)
nrow(y)

if (!is.null(cv)) {
  nrow(cv)
}
```

# 16. Important Notes About the Supplied Code

The manual describes the code as supplied rather than changing its statistical implementation.

A few implementation details should therefore be checked before using the code for a production analysis:

- `hoprr()` generates null genotypes with `rbinom(..., size = 2, prob = p_g)`, where `p_g` is the fitted genotype probability vector returned by the null model.
- Invalid permutation fits are discarded and additional permutations are generated until the requested number of valid permutations is obtained.
- `slfn()` relies on numerical Jacobians and can fail when the Jacobian is singular or poorly conditioned.
- The supplied `plafiu()` return list contains a duplicated name (`sd1`) for the two standard-error entries. This should be checked if the second standard error is intended to be named separately.
- The toy simulation and the analysis functions are contained in the same script; for large-scale real-data analysis, it is usually cleaner to separate simulation code from the reusable functions.

# 17. Reproducibility

The HoP-RR function accepts a random seed:

```{r, eval=FALSE}
seed = 123
```

Using the same input data, code, and seed makes the permutation analysis reproducible.

For a manuscript or software release, record:

- R version,
- code version/date,
- number of permutations,
- random seed,
- genotype coding,
- trait definitions,
- covariates,
- missing-data handling.

# 18. Quick Reference

| Function | Purpose |
|---|---|
| `hoprr()` | Main permutation-based HoP-RR test |
| `binomial_lrt()` | Computes the reverse-regression likelihood-ratio statistic |
| `lm_binom()` | Fits binomial logistic regression |
| `plafiu()` | Performs the two forward linear regressions |
| `marg_plei()` | Performs marginal pleiotropy tests |
| `numeric_gradient()` | Numerical gradient |
| `jacobian()` | Numerical Jacobian |
| `slfn()` | Newton-Raphson solver |

# 19. Minimal Example

Once the functions have been sourced, the simplest HoP-RR analysis is:

```{r, eval=FALSE}
source("hpl.R")

result <- hoprr(
  gt = my_genotype,
  y = cbind(my_trait1, my_trait2),
  cv = NULL,
  n_perm = 1000,
  seed = 123
)

result$p_value
result$decision
```

This is the recommended starting point for users who only want to apply the main HoP-RR function and do not need to modify the underlying numerical routines.

# 20. Conclusion

The supplied R script contains a complete implementation of the numerical routines, reverse-regression likelihood-ratio procedure, permutation-based HoP-RR test, forward-regression procedure, marginal test, and a toy simulation.

For routine analysis, users generally only need to prepare:

```text
genotype vector
two-trait matrix
optional covariate matrix
```

and call:

```r
hoprr(gt, y, cv)
```

The remaining functions support the internal estimation and testing procedure.
