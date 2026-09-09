## ==========================================================================
## Directed propensity-curvature contamination test (+ bias-corrected estimand)
##
## Null tested:  H0: theta_D = theta_Y = 0   (no consequential propensity
##               nonlinearity)  <=>  linear-covariate 2SLS keeps its
##               convex-weighted LATE interpretation.
##
## Key property we demonstrate: the test discriminates contaminated from
## clean designs EXACTLY where the Sanderson-Windmeijer conditional F cannot
## (F is ~enormous in both).
##
## Contamination moments (observable, no endogeneity, no weak-IV problem):
##   h(X) = E[Z|X] - L(Z|X)          nonlinear part of instrument propensity
##   theta_D = E[h(X) D],  theta_Y = E[h(X) Y]
## Bias-corrected estimand:  sum((Z-ehat)Y) / sum((Z-ehat)D)   (flexible ehat)
##   -- whose OWN first stage vanishes under collinearity (the no-free-lunch fact)
## ==========================================================================

set.seed(20260901)

## ---------- DGP ----------
gen <- function(n, prop = c("nonlinear", "linear"),
                kappa = 4, A_mu = 1.5, pi0 = 0.2, dp = 0.4) {
  prop <- match.arg(prop)
  X <- runif(n, -2, 2)
  g <- X^2 - 4/3
  e <- if (prop == "nonlinear") plogis(kappa * g)   # nonlinear -> h != 0
       else 0.5 + 0.1 * X                           # linear     -> h  = 0
  Z <- rbinom(n, 1, e)
  pi1 <- pi0 + dp; V <- runif(n)
  D  <- ifelse(Z == 1, as.integer(V <= pi1), as.integer(V <= pi0))
  tau <- 1 + 0.25 * X                # LATE(x) in [0.5,1.5], all positive
  Y   <- A_mu * g + tau * D + rnorm(n)
  list(X = X, Z = Z, D = D, Y = Y, tau = tau)
}

## ---------- one test on a dataset ----------
## Returns contamination estimates, Wald/p (pairs bootstrap), the naive 2SLS
## + conditional F, and the bias-corrected estimand + its own first-stage t.
contam_test <- function(d, degree = 3, B = 199) {
  X <- d$X; Z <- d$Z; D <- d$D; Y <- d$Y; n <- length(Z)
  P  <- poly(X, degree)             # cols: linear, quad, cubic (orthogonal)
  M  <- cbind(1, P)                 # full flexible design
  qi <- 2:degree                    # indices of curvature cols within P (quad, cubic)

  fit_h <- function(idx) {          # returns h_hat and ehat on rows idx
    Mi <- M[idx, , drop = FALSE]; Zi <- Z[idx]
    b  <- qr.solve(Mi, Zi)          # OLS of Z on [1, lin, quad, cubic]
    eh <- as.numeric(Mi %*% b)                       # flexible propensity
    hh <- as.numeric(P[idx, qi, drop = FALSE] %*% b[1 + qi]) # nonlinear part
    list(h = hh, e = eh)
  }

  f  <- fit_h(seq_len(n))
  thD <- mean(f$h * D); thY <- mean(f$h * Y)

  ## pairs bootstrap covariance of (thD, thY)
  bt <- matrix(NA, B, 2)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    fb  <- fit_h(idx)
    bt[b, ] <- c(mean(fb$h * D[idx]), mean(fb$h * Y[idx]))
  }
  Sig  <- cov(bt)
  th   <- c(thD, thY)
  Wald <- as.numeric(t(th) %*% solve(Sig, th))
  pval <- 1 - pchisq(Wald, df = 2)

  ## naive linear-covariate 2SLS + conditional (SW) first-stage F
  rz <- resid(lm(Z ~ X)); rd <- resid(lm(D ~ X)); ry <- resid(lm(Y ~ X))
  beta2sls <- sum(rz * ry) / sum(rz * rd)
  fs <- lm(D ~ Z + X); Fc <- coef(summary(fs))["Z", "t value"]^2

  ## bias-corrected estimand (flexible propensity residual) + its OWN strength
  zc   <- Z - f$e
  beta_corr <- sum(zc * Y) / sum(zc * D)
  t_corr    <- summary(lm(D ~ zc))$coefficients["zc", "t value"]  # corrected 1st-stage t

  c(thetaD = thD, thetaY = thY, Wald = Wald, pval = pval,
    beta2sls = beta2sls, Fcond = Fc, beta_corr = beta_corr, t_corr = t_corr)
}

## =================== (1) one large illustrative dataset ===================
cat("=== (1) Single contaminated dataset (n = 1e5) ===\n")
set.seed(1); d1 <- gen(1e5, "nonlinear", kappa = 4)
print(round(contam_test(d1, B = 199), 4))
cat("True LATE(x) support: [", round(min(d1$tau),2), ",", round(max(d1$tau),2), "]\n\n")

## =================== (2) size & power of the test ===================
mc <- function(prop, R = 300, n = 3000, B = 149) {
  rej <- Fc <- numeric(R)
  for (r in seq_len(R)) {
    out <- contam_test(gen(n, prop, kappa = 4), B = B)
    rej[r] <- out["pval"] < 0.05
    Fc[r]  <- out["Fcond"]
  }
  c(rej_rate_05 = mean(rej), mean_condF = mean(Fc))
}

cat("=== (2) Test properties (n = 3000, 5% level) ===\n")
set.seed(11); sz <- mc("linear")      # H0 true (propensity linear) -> SIZE
set.seed(22); pw <- mc("nonlinear")   # H0 false (curved propensity) -> POWER
tab <- rbind(`SIZE  (linear propensity, H0 true) ` = sz,
             `POWER (curved propensity, H0 false)` = pw)
print(round(tab, 3))
cat("\nBoth rows carry a HUGE conditional F, yet the test rejects only under\n",
    "genuine contamination -- the discrimination the conditional F cannot make.\n", sep = "")
