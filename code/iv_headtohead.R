## ==========================================================================
## Head-to-head: three tests of "linear-covariate 2SLS is uncontaminated",
## across three designs. Demonstrates (i) that first-stage-only tests cannot
## distinguish harmful from harmless propensity curvature (the impossibility
## result), and (ii) that the bias-targeting variant of our test is the only
## procedure correct in all three designs.
##
## Tests:
##   RESET_prop : F-test that the curvature basis has zero coefficients in the
##                propensity regression Z ~ [1, X, poly].  Uses only (Z, X).
##   Directed-S : joint Wald  H0: theta_D = theta_Y = 0   (sufficient condition;
##                strongly identified; cheap screening test).  chi^2_2.
##   Directed-X : bias-targeting moment  H0: E[h(X){Y - Lbar*D}] = 0, using the
##                bias-corrected Lbar-hat.  Necessary & sufficient for zero bias.
##                chi^2_1.
##
## Designs (n = 3000):
##   CONTAMINATED : nonlinear propensity + nonlinear outcome  -> real bias
##   HARMLESS     : nonlinear propensity + LINEAR outcome, const. effect -> NO bias
##   CLEAN        : linear propensity      + nonlinear outcome -> NO bias
## ==========================================================================

set.seed(20260902)

gen <- function(n, prop = c("nonlinear","linear"), outcome = c("nonlinear","linear"),
                kappa = 4, A_mu = 1.5, pi0 = 0.2, dp = 0.4) {
  prop <- match.arg(prop); outcome <- match.arg(outcome)
  X <- runif(n, -2, 2); g <- X^2 - 4/3
  e <- if (prop == "nonlinear") plogis(kappa * g) else pmin(pmax(0.5 + 0.1 * X, .02), .98)
  Z <- rbinom(n, 1, e)
  pi1 <- pi0 + dp; V <- runif(n)
  D <- ifelse(Z == 1, as.integer(V <= pi1), as.integer(V <= pi0))
  if (outcome == "nonlinear") { tau <- 1 + 0.25 * X; mu <- A_mu * g }  # heterog + curved level
  else                        { tau <- 1;            mu <- 0 }          # const effect + linear level
  Y <- mu + tau * D + rnorm(n)
  list(X = X, Z = Z, D = D, Y = Y)
}

## all three tests on one dataset (shared bootstrap for the two directed tests)
run_tests <- function(d, degree = 3, B = 99) {
  X <- d$X; Z <- d$Z; D <- d$D; Y <- d$Y; n <- length(Z)
  P <- poly(X, degree); M <- cbind(1, P); qi <- 2:degree

  ## RESET on the propensity (uses only Z, X): F-test of curvature block
  reset_p <- anova(lm(Z ~ X), lm(Z ~ P))[2, "Pr(>F)"]

  fit <- function(idx) {                 # returns h_hat, ehat, and Lbar on rows idx
    Mi <- M[idx, , drop = FALSE]; b <- qr.solve(Mi, Z[idx])
    eh <- as.numeric(Mi %*% b)
    hh <- as.numeric(P[idx, qi, drop = FALSE] %*% b[1 + qi])
    zc <- Z[idx] - eh
    Lbar <- sum(zc * Y[idx]) / sum(zc * D[idx])
    list(h = hh, Lbar = Lbar)
  }
  f0 <- fit(seq_len(n))
  thD <- mean(f0$h * D); thY <- mean(f0$h * Y)
  gX  <- mean(f0$h * (Y - f0$Lbar * D))            # bias-targeting moment

  bt <- matrix(NA, B, 3)                            # (thD, thY, gX)*
  for (k in seq_len(B)) {
    idx <- sample.int(n, n, TRUE); fk <- fit(idx)
    bt[k, ] <- c(mean(fk$h * D[idx]), mean(fk$h * Y[idx]),
                 mean(fk$h * (Y[idx] - fk$Lbar * D[idx])))
  }
  th   <- c(thD, thY)
  W_S  <- as.numeric(t(th) %*% solve(cov(bt[, 1:2]), th)); p_S <- 1 - pchisq(W_S, 2)
  W_X  <- gX^2 / var(bt[, 3]);                              p_X <- 1 - pchisq(W_X, 1)
  Fc   <- coef(summary(lm(D ~ Z + X)))["Z", "t value"]^2

  c(RESET = reset_p, DirS = p_S, DirX = p_X, Fcond = Fc)
}

mc <- function(prop, outcome, R = 300, n = 3000, B = 99) {
  M <- t(replicate(R, run_tests(gen(n, prop, outcome), B = B)))
  c(RESET = mean(M[, "RESET"] < .05),
    Directed_S = mean(M[, "DirS"] < .05),
    Directed_X = mean(M[, "DirX"] < .05),
    mean_condF = mean(M[, "Fcond"]))
}

cat("Running head-to-head (3 designs x 3 tests)...\n")
set.seed(1);  con <- mc("nonlinear", "nonlinear")   # bias != 0
set.seed(2);  har <- mc("nonlinear", "linear")      # bias  = 0 (harmless curvature)
set.seed(3);  cln <- mc("linear",    "nonlinear")   # bias  = 0 (clean)

tab <- rbind(`CONTAMINATED (bias != 0)` = con,
             `HARMLESS     (bias  = 0)` = har,
             `CLEAN        (bias  = 0)` = cln)
cat("\n=== Rejection rates (5% level, R=300, n=3000) ===\n")
print(round(tab, 3))
cat("\nReading:\n",
    " CONTAMINATED : all three should reject (power).\n",
    " HARMLESS     : RESET and Directed-S FALSE-ALARM (they flag curvature, not bias);\n",
    "                Directed-X stays at size -- the only test that certifies a fine design.\n",
    " CLEAN        : all three at size.\n",
    " mean_condF is huge in EVERY row -> the first stage cannot tell these apart.\n", sep = "")

## ------------------------------- figure -------------------------------
png("iv_headtohead.png", width = 1450, height = 620, res = 150)
op <- par(mar = c(4.5, 4.6, 3, 1))
M <- t(tab[, 1:3])                                  # tests x designs
cols <- c("#b0b0b0", "#7570b3", "#d95f02")
bp <- barplot(M, beside = TRUE, col = cols, ylim = c(0, 1.08),
              ylab = "rejection rate (5% level)",
              names.arg = c("CONTAMINATED\n(bias ≠ 0)", "HARMLESS\n(bias = 0)", "CLEAN\n(bias = 0)"),
              main = "Three tests across three designs")
abline(h = 0.05, lty = 2, col = "gray40")
legend("topright", bty = "n",
       legend = c("RESET (propensity only; ignores Y)",
                  "Directed-S (sufficient: theta=0)",
                  "Directed-X (bias-targeting)"),
       fill = cols, cex = 0.8)
text(mean(bp[, 2]), 0.9, "RESET & Directed-S\nfalse-alarm here", cex = 0.7, col = "gray25")
par(op); dev.off()
cat("\nSaved figure -> iv_headtohead.png\n")
