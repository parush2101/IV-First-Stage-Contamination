## ==========================================================================
## Monte Carlo test of the "collinearity contamination" proposition
##
## Claim: In just-identified 2SLS with covariates entered LINEARLY
## (non-saturated), as the instrument Z becomes more determined by X
## (instrument-covariate collinearity), the 2SLS estimand departs
## arbitrarily from every conditional LATE(x) -- WHILE the Sanderson-Windmeijer
## conditional first-stage F stays large. I.e. the conditional F fails to warn.
##
## Structural DGP (conditional-on-X IV validity + conditional monotonicity):
##   Y = mu(X) + tau(X)*D + noise      tau(X)=LATE(X) heterogeneous, all > 0
##   D = D(Z),  monotone potential treatments (compliers/always/never takers)
##   e(X)=P(Z=1|X) nonlinear in X;  kappa = collinearity knob (Var(Z|X)->0)
##
## We report, across kappa:
##   beta2sls      : linear-covariate 2SLS estimand (finite-sample, large N)
##   Fcond         : conditional (SW) first-stage F, HC1-robust
##   late_convex   : the convex-weighted conditional LATE it SHOULD be near
##   [late_min,max]: support of the true conditional effects
##   contam_limit  : theoretical limit Cov(h,mY)/Cov(h,mD)  (h = nonlinear part of e)
## ==========================================================================

set.seed(20260901)

## HC1-robust t-statistic for one named regressor in an lm fit
robust_t <- function(fit, coef_name) {
  X  <- model.matrix(fit)
  u  <- resid(fit)
  n  <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X))
  meat  <- t(X) %*% (X * (u^2))
  V     <- bread %*% meat %*% bread * (n / (n - k))   # HC1
  b     <- coef(fit)[coef_name]
  se    <- sqrt(V[coef_name, coef_name])
  as.numeric(b / se)
}

simulate <- function(N, kappa,
                     A_mu = 1.5,    # strength of OUTCOME nonlinearity in X
                     pi0  = 0.2,    # P(always-taker | x)
                     dp   = 0.4) {  # conditional first stage p(x)=P(complier|x)
  X <- runif(N, -2, 2)
  g <- X^2 - (4 / 3)               # E[X^2]=4/3 on U(-2,2): centered nonlinear shape
  e <- plogis(kappa * g)           # instrument propensity, nonlinear in X
  Z <- rbinom(N, 1, e)

  ## monotone potential treatments (pi1 >= pi0 => no defiers)
  pi1 <- pi0 + dp
  V   <- runif(N)
  D0  <- as.integer(V <= pi0)
  D1  <- as.integer(V <= pi1)
  D   <- ifelse(Z == 1, D1, D0)

  ## heterogeneous, strictly positive effect: LATE(x) in [0.5, 1.5]
  tau <- 1 + 0.25 * X
  mu  <- A_mu * g                  # nonlinear covariate level (shares shape with e)
  Y   <- mu + tau * D + rnorm(N, 0, 1)

  ## linear-covariate (NON-saturated) 2SLS via Frisch-Waugh-Lovell
  rz <- resid(lm(Z ~ X))
  rd <- resid(lm(D ~ X))
  ry <- resid(lm(Y ~ X))
  beta2sls <- sum(rz * ry) / sum(rz * rd)

  ## conditional (Sanderson-Windmeijer) first-stage F = robust t^2 on Z
  fs    <- lm(D ~ Z + X)
  Fcond <- robust_t(fs, "Z")^2

  ## target the estimand SHOULD approximate: convex-weighted conditional LATE
  w           <- e * (1 - e) * dp          # signal weights e(1-e)p(x)
  late_convex <- sum(w * tau) / sum(w)

  ## theoretical contamination limit Cov(h,mY)/Cov(h,mD)
  h  <- resid(lm(e ~ X))                    # nonlinear part of the propensity
  mD <- pi0 + e * dp                        # E[D|X]
  mY <- mu + tau * mD                       # E[Y|X]
  contam_limit <- cov(h, mY) / cov(h, mD)

  c(kappa = kappa, beta2sls = beta2sls, Fcond = Fcond,
    late_convex = late_convex, late_min = min(tau), late_max = max(tau),
    contam_limit = contam_limit, ownvar_ZgivenX = mean(e * (1 - e)))
}

N      <- 200000
kappas <- c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
res    <- as.data.frame(t(sapply(kappas, function(k) simulate(N, k))))

cat("\n=== IV collinearity-contamination Monte Carlo (N =", N, "per row) ===\n\n")
print(round(res, 3), row.names = FALSE)

## ---------------------------- figure ----------------------------
png("iv_collinearity_fig.png", width = 1500, height = 640, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3, 1.2))

## Panel 1: the estimand vs the truth
plot(res$kappa, res$beta2sls, type = "b", pch = 19, log = "x",
     xlab = expression(paste("collinearity  ", kappa, "   (", Var(Z*"|"*X), " ", symbol("\256"), " 0)")),
     ylab = "estimand", main = "2SLS estimand leaves the LATE hull",
     ylim = range(0, res$beta2sls, res$late_max) * c(1, 1.05))
polygon(c(res$kappa, rev(res$kappa)),
        c(res$late_min, rev(res$late_max)),
        col = rgb(0.2, 0.6, 0.2, 0.18), border = NA)
lines(res$kappa, res$late_convex, col = "forestgreen", lwd = 2, lty = 2)
lines(res$kappa, res$contam_limit, col = "firebrick", lwd = 1.5, lty = 3)
lines(res$kappa, res$beta2sls, type = "b", pch = 19, lwd = 2)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("2SLS estimand", "true convex-weighted LATE",
                  "band of all LATE(x)", "contamination limit"),
       col = c("black", "forestgreen", rgb(0.2,0.6,0.2,0.5), "firebrick"),
       lwd = c(2, 2, 8, 1.5), lty = c(1, 2, 1, 3), pch = c(19, NA, NA, NA))

## Panel 2: the diagnostic that fails to warn
plot(res$kappa, res$Fcond, type = "b", pch = 19, log = "xy", col = "firebrick",
     xlab = expression(paste("collinearity  ", kappa)),
     ylab = "conditional first-stage F (log scale)",
     main = "...while the conditional F stays enormous")
abline(h = 10, lty = 2, col = "gray40")
text(min(res$kappa), 13, "weak-IV threshold F = 10", pos = 4, cex = 0.8, col = "gray40")
par(op); dev.off()
cat("\nSaved figure -> iv_collinearity_fig.png\n")

cat("\nRead the table top-to-bottom as collinearity (kappa) rises:\n")
cat(" * ownvar_ZgivenX -> 0   : Z becomes ~deterministic in X (the collinearity limit)\n")
cat(" * Fcond stays HUGE       : the conditional first stage looks rock-solid throughout\n")
cat(" * beta2sls leaves [late_min, late_max]: estimand exits the hull of ALL true effects\n")
cat(" * beta2sls -> contam_limit: it converges to the pure-contamination ratio, not a LATE\n")
