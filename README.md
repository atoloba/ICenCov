The final version of the package will be available by November 2025

# ICenCov: Regression Models with Interval-Censored Covariates 

Function **icglm** fits Generalized Linear Models with an interval-censored covariate.
It uses the **GELc** estimator: a semiparametric, likelihood-based approach 
built on an augmented version of Turnbull's nonparametric estimator for 
interval-censored data. 

## Installation

**ICenCov** can be installed from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("atoloba/ICenCov")
```

## Brief Examples

### Logistic regression

We use the *carotenoids* dataset, included in this package, to illustrate how to fit a logistic regression model 
with an interval-censored covariate. The variables included in the model are:

- *diabetes status* as response variable
- *sex* and *bmi* as adjustment covariates
- *caro* as covariate of interest, interval-censored in [*zl*, *zr*)


``` r
# load the dataset
data(carotenoids)
dcaro <- subset(carotenoids, !is.na(diabetes))

# fit the logistic regression
fit1 <- icglm(diabetes ~ sex + bmi + ic(zl, zr, "caro"), family = binomial,
              data = dcaro, Lin = TRUE,
              Rin  = with(dcaro, zl == zr)
             )

# print the summary
summary(fit1)

# Odds ratio (OR) and 95% CI for 'caro':
i = 4
OR_caro  <- exp(coef(fit1)[i])
CI_caro  <- exp(coef(fit1)[i] + c(-1,1)*qnorm(0.975)*sqrt(vcov(fit1)[i,i]))
cat(sprintf("OR = %.2f (%.2f, %.2f)\n", OR_caro, CI_caro[1], CI_caro[2]))
```

### Gamma regression model

Using the ACTG359 dataset included in the package, we show how to fit a generalized linear model
under the *Gamma(link="log")* family. The expression of the model is
```math
E[{\rm RNA}] = \exp\{ \alpha+\beta {\rm age} + \gamma {\rm waitime} \},
```
and *waitime* is interval-censored in (*zl*, *zr*].

<br>

``` r
# load the dataset
data(actg359)

# fit the gamma regression model specifying the link function 'log'
fit2 <- icglm(RNA ~ age + ic(zl, zr, "waitime"), family = Gamma("log"),
              data = actg359, Lin = FALSE, Rin = TRUE)

# print the summary
summary(fit2)

# Multiplicative effect of a 3-week increase in 'waitime', and 95% CI
i = 3
mult_wt   <- exp(3 * coef(fit2)[i]) - 1
CI_wt <- exp(3 * coef(fit2)[i] + 3*c(-1,1)*qnorm(0.975)*sqrt(vcov(fit2)[i,i])) - 1
cat(sprintf("OR = %.1f (%.1f, %.1f)\n", mult_wt*100, CI_wt[1]*100, CI_wt[2]*100))
```


## Citation

```bibtex
@misc{ICenCov,
      title  = {{IC}en{C}ov: {R}egression {M}odels with {I}nterval-{C}ensored {C}ovariates}, 
      author = {Andrea Toloba and Klaus Langohr and Guadalupe {Gómez Melis}},
      year   = {2025},
      url    = {https://github.com/atoloba/ICenCov}
}
```
