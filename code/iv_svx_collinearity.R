## ==========================================================================
## Task 6: Directed-S vs Directed-X as collinearity rises (a CONTAMINATED design
## with GENUINE bias throughout). Substantiates the "report both" recommendation:
##   - Directed-S (sufficient screen, theta=0) never forms an IV ratio, so it
##     stays strongly identified and keeps full power as Var(Z|X) -> 0.
##   - Directed-X (exact, bias-targeting) plugs in the corrected Lbar-hat, which
##     is weakly identified as S_D -> 0 (Prop 6 / Cor). Its studentization inflates
##     and it LOSES power under severe collinearity, even though the bias grows.
## Same DGP and tests as iv_headtohead.R; contaminated design (nonlinear prop +
## nonlinear outcome), swept over the collinearity knob kappa.
## ==========================================================================
set.seed(20260902)

gen <- function(n, kappa, A_mu = 1.5, pi0 = 0.2, dp = 0.4) {
  X <- runif(n, -2, 2); g <- X^2 - 4/3
  e <- plogis(kappa * g)
  Z <- rbinom(n, 1, e)
  pi1 <- pi0 + dp; V <- runif(n)
  D <- ifelse(Z == 1, as.integer(V <= pi1), as.integer(V <= pi0))
  tau <- 1 + 0.25 * X; mu <- A_mu * g
  Y <- mu + tau * D + rnorm(n)
  list(X = X, Z = Z, D = D, Y = Y)
}

run_tests <- function(d, degree = 3, B = 99) {
  X <- d$X; Z <- d$Z; D <- d$D; Y <- d$Y; n <- length(Z)
  P <- poly(X, degree); M <- cbind(1, P); qi <- 2:degree
  reset_p <- anova(lm(Z ~ X), lm(Z ~ P))[2, "Pr(>F)"]
  fit <- function(idx) {
    Mi <- M[idx, , drop = FALSE]; b <- qr.solve(Mi, Z[idx])
    eh <- as.numeric(Mi %*% b)
    hh <- as.numeric(P[idx, qi, drop = FALSE] %*% b[1 + qi])
    zc <- Z[idx] - eh
    list(h = hh, Lbar = sum(zc * Y[idx]) / sum(zc * D[idx]))
  }
  f0 <- fit(seq_len(n))
  thD <- mean(f0$h * D); thY <- mean(f0$h * Y)
  gX  <- mean(f0$h * (Y - f0$Lbar * D))
  bt <- matrix(NA, B, 3)
  for (k in seq_len(B)) {
    idx <- sample.int(n, n, TRUE); fk <- fit(idx)
    bt[k, ] <- c(mean(fk$h * D[idx]), mean(fk$h * Y[idx]),
                 mean(fk$h * (Y[idx] - fk$Lbar * D[idx])))
  }
  th <- c(thD, thY)
  W_S <- as.numeric(t(th) %*% solve(cov(bt[, 1:2]), th)); p_S <- 1 - pchisq(W_S, 2)
  W_X <- gX^2 / var(bt[, 3]);                             p_X <- 1 - pchisq(W_X, 1)
  Fc  <- coef(summary(lm(D ~ Z + X)))["Z", "t value"]^2
  c(RESET = reset_p, DirS = p_S, DirX = p_X, Fcond = Fc)
}

## population bias and contamination share at large N (for context columns)
popinfo <- function(kappa, N = 200000, A_mu = 1.5, pi0 = 0.2, dp = 0.4) {
  X <- runif(N, -2, 2); g <- X^2 - 4/3; e <- plogis(kappa * g); Z <- rbinom(N, 1, e)
  pi1 <- pi0 + dp; V <- runif(N)
  D <- ifelse(Z == 1, as.integer(V <= pi1), as.integer(V <= pi0))
  tau <- 1 + 0.25 * X; mu <- A_mu * g; Y <- mu + tau * D + rnorm(N)
  rz <- resid(lm(Z ~ X)); rd <- resid(lm(D ~ X)); ry <- resid(lm(Y ~ X))
  beta <- sum(rz * ry) / sum(rz * rd)
  w <- e * (1 - e) * dp; Lconv <- sum(w * tau) / sum(w)
  h <- resid(lm(e ~ X)); thD <- mean(h * D); rho <- thD / mean(rz * D)
  c(bias = beta - Lconv, rho = rho, ownvar = mean(e * (1 - e)))
}

mc <- function(kappa, R = 300, n = 3000, B = 99) {
  MM <- t(replicate(R, run_tests(gen(n, kappa), B = B)))
  pil <- popinfo(kappa)
  c(kappa = kappa, rho = pil["rho"], bias = pil["bias"],
    RESET = mean(MM[, "RESET"] < .05),
    Directed_S = mean(MM[, "DirS"] < .05),
    Directed_X = mean(MM[, "DirX"] < .05),
    mean_condF = mean(MM[, "Fcond"]))
}

kappas <- c(0.5, 1, 2, 4, 8, 16)
cat("Running Directed-S vs Directed-X sweep over collinearity (contaminated design)...\n")
tab <- as.data.frame(t(sapply(seq_along(kappas), function(i){ set.seed(100 + i); mc(kappas[i]) })))
names(tab) <- c("kappa","rho","bias","RESET","Directed_S","Directed_X","mean_condF")
cat("\n=== Rejection rates (5% level, R=300, n=3000), contaminated design ===\n")
print(round(tab, 3), row.names = FALSE)
cat("\nReading: bias GROWS with kappa, yet Directed_X power FALLS at high kappa\n",
    "(corrected Lbar weakly identified), while Directed_S stays ~1 throughout.\n", sep="")
