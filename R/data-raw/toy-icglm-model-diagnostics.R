## data-raw/toy-icglm-model-diagnostics.R
## ------------------------------------------------------------------
## Precompute toy datasets + fitted models for fast @examples
## ------------------------------------------------------------------
devtools::load_all()
library(dplyr)

## ================================================================
## 1) Toy example: residual normality check
## ================================================================
set.seed(123)
n <- 300
alpha <- 1
gamma <- 0.5
phi <- 0.5

z_true <- runif(n, 0, 5)
mu <- exp(alpha + gamma * z_true)
y <- rgamma(n, shape = 1/phi, scale = mu * phi)

mgap <- 0.5
sdgap <- 0.75 * mgap

cens <- censoring(
  z_true,
  rini = function(nn) runif(nn, 0, mgap),
  rfun = function(nn) rnorm(nn, mgap, sdgap)
)
if (!is.null(cens$warning)) message(cens$warning)

dat <- data.frame(
  y = y,
  zl = cens$zl,
  zr = cens$zr,
  z_true = z_true
)

fit_gam <- icglm(
  y ~ ic(zl, zr),
  Gamma(link = "log"),
  data = dat,
  tolik = 1e-6
)

fit_ig <- icglm(
  y ~ ic(zl, zr),
  inverse.gaussian(link = "log"),
  data = dat,
  tolik = 1e-6
)

toy_resnorm <- list(dat = dat, fit_gam = fit_gam, fit_ig = fit_ig)
usethis::use_data(toy_resnorm, overwrite = TRUE)


## ================================================================
## 2) Toy example: linearity of covariates
## ================================================================
set.seed(123)
n <- 300

z_true <- runif(n, 0, 50)  # higher range because we will take log(z)
x <- rnorm(n, 0, 0.5)

alpha <- 0
beta  <- 0.5   # U-shape -> suggests x^2
gamma <- 1     # concave -> suggests log(z)
phi   <- 0.01

mu <- exp(alpha + beta * x^2 + gamma * log(z_true))
y <- rgamma(n, shape = 1/phi, scale = mu * phi)

mgap <- 5
sdgap <- 0.75 * mgap

cens <- censoring(
  z_true,
  rini = function(nn) runif(nn, 0, mgap),
  rfun = function(nn) rnorm(nn, mgap, sdgap)
)
if (!is.null(cens$warning)) message(cens$warning)

dat <- data.frame(
  y = y,
  zl = cens$zl,
  zr = cens$zr,
  z_true = z_true,
  x = x
)

## “Naive” model: linear x and z
fit_linear <- icglm(
  y ~ x + ic(zl, zr),
  Gamma(link = "log"),
  data = dat,
  tolik = 1e-6
)

## Transformed covariates: x^2 and log z
dat <- dat %>% mutate(
  x2  = x^2,
  zl2 = log(pmax(zl, 0.03)),
  zr2 = log(zr)
)

fit_trans <- icglm(
  y ~ x2 + ic(zl2, zr2),
  Gamma(link = "log"),
  data = dat,
  tolik = 1e-6
)

toy_linearity <- list(dat = dat, fit_linear = fit_linear, fit_trans = fit_trans)
usethis::use_data(toy_linearity, overwrite = TRUE)


## ================================================================
## 3) Toy example: link function check (identity vs log)
## ================================================================
set.seed(123)
n <- 300
alpha <- 0
gamma <- 0.5
phi <- 0.1

z_true <- runif(n, 0, 5)
mu <- exp(alpha + gamma * z_true)
y <- rgamma(n, shape = 1/phi, scale = mu * phi)

mgap <- 0.5
sdgap <- 0.75 * mgap

cens <- censoring(
  z_true,
  rini = function(nn) runif(nn, 0, mgap),
  rfun = function(nn) rnorm(nn, mgap, sdgap)
)
if (!is.null(cens$warning)) message(cens$warning)

dat <- data.frame(
  y = y,
  zl = cens$zl,
  zr = cens$zr,
  z_true = z_true
)

fit_id <- icglm(
  y ~ ic(zl, zr),
  Gamma(link = "identity"),
  data = dat,
  tolik = 1e-6
)

fit_log <- icglm(
  y ~ ic(zl, zr),
  Gamma(link = "log"),
  data = dat,
  tolik = 1e-6
)

toy_link <- list(dat = dat, fit_id = fit_id, fit_log = fit_log)
usethis::use_data(toy_link, overwrite = TRUE)
