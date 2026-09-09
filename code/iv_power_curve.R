## ==========================================================================
## Power curve for the directed propensity-curvature contamination test.
##
## Contamination knob b: instrument propensity is a LINEAR-probability model
##   e(X) = 0.5 + a*X + b*(X^2 - 4/3),  clipped to (0.02, 0.98)
##   b = 0  => propensity EXACTLY linear in X => h = 0 => H0 true (no bias)
##   b > 0  => curvature => contamination grows monotonically
##
## We plot rejection rate vs b for three sample sizes, and (2nd panel) the
## mean conditional (Sanderson-Windmeijer) F vs b -- which stays far above the
## weak-IV threshold across the ENTIRE curve, i.e. gives a green light where
## the test correctly escalates from size to power.
## ==========================================================================

set.seed(20260901)

gen_lpm <- function(n, b, a = 0.08, A_mu = 1.5, pi0 = 0.2, dp = 0.4) {
  X <- runif(n, -2, 2); g <- X^2 - 4/3
  e <- pmin(pmax(0.5 + a * X + b * g, 0.02), 0.98)
  Z <- rbinom(n, 1, e)
  pi1 <- pi0 + dp; V <- runif(n)
  D <- ifelse(Z == 1, as.integer(V <= pi1), as.integer(V <= pi0))
  tau <- 1 + 0.25 * X
  Y <- A_mu * g + tau * D + rnorm(n)
  list(X = X, Z = Z, D = D, Y = Y, tau = tau)
}

## lean test: bootstrap p-value + conditional F
test_pval <- function(d, degree = 3, B = 99) {
  X <- d$X; Z <- d$Z; D <- d$D; Y <- d$Y; n <- length(Z)
  P <- poly(X, degree); M <- cbind(1, P); qi <- 2:degree
  hfit <- function(idx) {
    Mi <- M[idx, , drop = FALSE]
    bb <- qr.solve(Mi, Z[idx])
    as.numeric(P[idx, qi, drop = FALSE] %*% bb[1 + qi])
  }
  h  <- hfit(seq_len(n)); th <- c(mean(h * D), mean(h * Y))
  bt <- matrix(NA, B, 2)
  for (k in seq_len(B)) {
    idx <- sample.int(n, n, TRUE); hk <- hfit(idx)
    bt[k, ] <- c(mean(hk * D[idx]), mean(hk * Y[idx]))
  }
  Wald <- as.numeric(t(th) %*% solve(cov(bt), th))
  p    <- 1 - pchisq(Wald, 2)
  Fc   <- coef(summary(lm(D ~ Z + X)))["Z", "t value"]^2
  c(p = p, Fc = Fc)
}

mc_cell <- function(n, b, R = 200, B = 99) {
  pr <- Fc <- numeric(R)
  for (r in seq_len(R)) { o <- test_pval(gen_lpm(n, b), B = B); pr[r] <- o["p"]; Fc[r] <- o["Fc"] }
  c(rej = mean(pr < 0.05), F = mean(Fc))
}

## population bias map: b -> |beta2sls - convex LATE|  (one big-n draw)
bias_at <- function(b, n = 2e5) {
  d <- gen_lpm(n, b)
  rz <- resid(lm(d$Z ~ d$X)); rd <- resid(lm(d$D ~ d$X)); ry <- resid(lm(d$Y ~ d$X))
  beta <- sum(rz * ry) / sum(rz * rd)
  e  <- pmin(pmax(0.5 + 0.08 * d$X + b * (d$X^2 - 4/3), 0.02), 0.98)
  w  <- e * (1 - e); Lbar <- sum(w * d$tau) / sum(w)
  beta - Lbar
}

bs   <- seq(0, 0.12, by = 0.02)
ns   <- c(1000, 3000, 9000)
cat("Running power grid (", length(bs) * length(ns), "cells)...\n", sep = "")

rej <- Fm <- matrix(NA, length(bs), length(ns), dimnames = list(paste0("b=", bs), paste0("n=", ns)))
for (j in seq_along(ns)) for (i in seq_along(bs)) {
  set.seed(100 * j + i)
  o <- mc_cell(ns[j], bs[i]); rej[i, j] <- o["rej"]; Fm[i, j] <- o["F"]
}
set.seed(7); bias <- sapply(bs, bias_at)

cat("\n=== Rejection rate (5% level) by contamination knob b and n ===\n")
print(round(rej, 3))
cat("\n=== Population bias  beta2sls - convex-LATE  at each b ===\n")
print(round(setNames(bias, paste0("b=", bs)), 3))
cat("\n=== Mean conditional F (n=9000 column) ===\n")
print(round(Fm[, 3], 0))

## ------------------------------- figure -------------------------------
png("iv_power_curve.png", width = 1500, height = 640, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3, 1.2))
cols <- c("#1b9e77", "#7570b3", "#d95f02")

plot(bs, rej[, 1], type = "b", pch = 19, col = cols[1], ylim = c(0, 1),
     xlab = "propensity curvature  b  (contamination magnitude)",
     ylab = "rejection rate  (5% level)",
     main = "Power of the contamination test")
for (j in 2:3) lines(bs, rej[, j], type = "b", pch = 19, col = cols[j], lwd = 2)
lines(bs, rej[, 1], type = "b", pch = 19, col = cols[1], lwd = 2)
abline(h = 0.05, lty = 2, col = "gray50")
axis(3, at = bs, labels = round(bias, 1), cex.axis = 0.7)
mtext("population bias  (β − convex-LATE)", side = 3, line = 2, cex = 0.8)
legend("right", bty = "n", legend = paste0("n = ", ns), col = cols, lwd = 2, pch = 19)
text(0, 0.11, "size at b=0", pos = 4, cex = 0.75, col = "gray40")

plot(bs, Fm[, 3], type = "b", pch = 19, log = "y", col = "firebrick", lwd = 2,
     xlab = "propensity curvature  b", ylab = "mean conditional F (log scale)",
     main = "...the conditional F says 'fine' throughout", ylim = c(5, max(Fm)))
abline(h = 10, lty = 2, col = "gray40")
text(0, 12, "weak-IV threshold F = 10", pos = 4, cex = 0.8, col = "gray40")
par(op); dev.off()
cat("\nSaved figure -> iv_power_curve.png\n")
