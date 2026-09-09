## ==========================================================================
## The contamination share of the first stage.
##
##   rho  =  theta_D / E[Ztilde D]  =  theta_D / (S_D + theta_D)
##        =  fraction of the observed first-stage covariance (hence of the
##           apparent strength) that comes from propensity-curvature
##           contamination rather than genuine complier variation.
##
## In the nonlinear-propensity design this equals E[h^2]/Var(Ztilde), the share
## of the residualized-instrument variance that is curvature.  We show the
## plug-in estimate is consistent, its percentile-bootstrap CI covers, and rho
## rises toward 1 exactly as the 2SLS bias explodes -- a readable danger gauge
## to report next to the conditional F.
## ==========================================================================

set.seed(20260902)

gen <- function(n, kappa, A_mu = 1.5, pi0 = 0.2, dp = 0.4) {
  X <- runif(n, -2, 2); g <- X^2 - 4/3
  e <- plogis(kappa * g); Z <- rbinom(n, 1, e)
  pi1 <- pi0 + dp; V <- runif(n)
  D <- ifelse(Z == 1, as.integer(V <= pi1), as.integer(V <= pi0))
  tau <- 1 + 0.25 * X
  Y <- A_mu * g + tau * D + rnorm(n)
  list(X = X, Z = Z, D = D, Y = Y, tau = tau)
}

## plug-in contamination share on rows idx (design matrix M, basis P precomputable)
share_hat <- function(X, Z, D, degree = 3) {
  P <- poly(X, degree); M <- cbind(1, P); qi <- 2:degree
  b  <- qr.solve(M, Z)
  hh <- as.numeric(P[, qi, drop = FALSE] %*% b[1 + qi])
  rz <- resid(lm(Z ~ X))
  C   <- mean(rz * D)          # E[Ztilde D]  (the SW first-stage numerator)
  thD <- mean(hh * D)          # contamination part
  thD / C
}

boot_ci <- function(d, B = 149, degree = 3) {
  n <- length(d$Z)
  bs <- numeric(B)
  for (k in seq_len(B)) {
    idx <- sample.int(n, n, TRUE)
    bs[k] <- share_hat(d$X[idx], d$Z[idx], d$D[idx], degree)
  }
  quantile(bs, c(.025, .975), names = FALSE)
}

kappas <- c(0.25, 0.5, 1, 2, 4, 8, 16)

## ---- (A) large-n truth: true share and true bias per kappa ----
cat("Computing large-n truth...\n")
truth <- t(sapply(kappas, function(k) {
  d <- gen(2e5, k)
  rho_true <- share_hat(d$X, d$Z, d$D)
  rz <- resid(lm(d$Z ~ d$X)); rd <- resid(lm(d$D ~ d$X)); ry <- resid(lm(d$Y ~ d$X))
  beta <- sum(rz * ry) / sum(rz * rd)
  e <- plogis(k * (d$X^2 - 4/3)); w <- e * (1 - e); Lbar <- sum(w * d$tau) / sum(w)
  c(rho_true = rho_true, bias = beta - Lbar)
}))

## ---- (B) moderate-n: mean estimate + CI coverage of the true share ----
cat("Running moderate-n Monte Carlo (n=3000)...\n")
R <- 250; n <- 3000
modres <- t(sapply(seq_along(kappas), function(j) {
  k <- kappas[j]; rho0 <- truth[j, "rho_true"]
  est <- cover <- numeric(R)
  for (r in seq_len(R)) {
    d <- gen(n, k)
    est[r] <- share_hat(d$X, d$Z, d$D)
    ci <- boot_ci(d, B = 149)
    cover[r] <- (rho0 >= ci[1] & rho0 <= ci[2])
  }
  c(mean_est = mean(est), coverage95 = mean(cover))
}))

tab <- data.frame(kappa = kappas,
                  rho_true = round(truth[, "rho_true"], 3),
                  mean_est_n3000 = round(modres[, "mean_est"], 3),
                  CI95_coverage = round(modres[, "coverage95"], 3),
                  bias = round(truth[, "bias"], 3))
cat("\n=== Contamination share: truth, estimate, coverage, and the induced bias ===\n")
print(tab, row.names = FALSE)

## ------------------------------- figure -------------------------------
## one moderate-n dataset per kappa for a CI band illustration
set.seed(99)
band <- t(sapply(kappas, function(k) { d <- gen(3000, k); c(est = share_hat(d$X,d$Z,d$D), boot_ci(d, 199)) }))

png("iv_contamination_share.png", width = 1450, height = 620, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3, 1.1))

plot(kappas, truth[, "rho_true"], type = "l", lwd = 2, col = "firebrick", log = "x",
     ylim = c(0, 1), xlab = expression(paste("collinearity  ", kappa)),
     ylab = expression(paste("contamination share  ", varrho)),
     main = "Contamination share (truth, estimate, 95% CI)")
polygon(c(kappas, rev(kappas)), c(band[,2], rev(band[,3])),
        col = rgb(.2,.2,.6,.15), border = NA)
points(kappas, band[,1], pch = 19); lines(kappas, band[,1], lty = 3)
legend("bottomright", bty = "n", cex = .8,
       legend = c("true share", "plug-in estimate (n=3000)", "95% bootstrap CI"),
       col = c("firebrick","black",rgb(.2,.2,.6,.5)), lwd = c(2,1,8), lty = c(1,3,1), pch = c(NA,19,NA))

plot(truth[, "rho_true"], truth[, "bias"], type = "b", pch = 19, lwd = 2,
     xlab = expression(paste("contamination share  ", varrho)),
     ylab = expression(paste(beta[2*SLS], " - convex-LATE")),
     main = "Share is a readable proxy for the bias")
par(op); dev.off()
cat("\nSaved figure -> iv_contamination_share.png\n")
