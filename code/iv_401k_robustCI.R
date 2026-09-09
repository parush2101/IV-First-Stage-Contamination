# Identification-robust inference for the contamination-corrected estimand
# in the 401(k) linear-income design (closes the loop on Proposition 5 /
# recommendation (3)).  Corrected estimand Lbar is a just-identified IV with
# instrument Zc = Z - ehat(X); Anderson-Rubin is exact and identification-robust
# for just-identified IV, so we invert it for a CI and compare to the pairs
# bootstrap.  We also report the corrected estimator's OWN first-stage F (a
# sample analogue of S_D), which certifies strong identification in THIS design.
set.seed(20260902)
d <- foreign::read.dta("401ksubs.dta")

Y <- d$nettfa; D <- d$p401k; Z <- d$e401k
Xl <- with(d, cbind(inc, age, marr, male, fsize))                 # linear controls
Bt <- resid(lm(with(d, cbind(inc^2, inc^3, age^2, inc*age)) ~ Xl))# curvature basis
n  <- length(Z)

## corrected instrument Zc = Z - ehat(X), ehat flexible (linear controls + curvature)
ehat <- fitted(lm(Z ~ Xl + Bt))
Zc   <- Z - ehat
Lbar <- sum(Zc*Y) / sum(Zc*D)

## corrected estimator's OWN first stage: regress D on Zc (+controls), robust F
fs   <- lm(D ~ Zc + Xl)
X1   <- model.matrix(fs); u1 <- resid(fs); k1 <- ncol(X1)
br   <- solve(crossprod(X1)); V1 <- br %*% (t(X1)%*%(X1*u1^2)) %*% br * (n/(n-k1))
F_corr <- (coef(fs)["Zc"]/sqrt(V1["Zc","Zc"]))^2

## Anderson-Rubin confidence set by grid inversion.
## Moment m_i(b) = Zc_i (Y_i - b D_i);  AR(b) = (sum m)^2 / sum(m^2)  ~ chi^2_1.
AR <- function(b){ m <- Zc*(Y - b*D); (sum(m))^2 / sum(m^2) }
grid <- seq(-40, 60, by = 0.01)
arv  <- sapply(grid, AR)
inset <- grid[arv <= qchisq(0.95, 1)]
ar_lo <- min(inset); ar_hi <- max(inset)

## Wald (delta-method) SE for reference
mm <- Zc*(Y - Lbar*D); g <- sum(Zc*D)
se_wald <- sqrt(sum(mm^2))/abs(g)     # robust SE of just-identified IV
w_lo <- Lbar - 1.96*se_wald; w_hi <- Lbar + 1.96*se_wald

cat(sprintf("corrected Lbar         = %.2f\n", Lbar))
cat(sprintf("corrected own first-stage F (D ~ Zc) = %.0f\n", F_corr))
cat(sprintf("Wald 95%% CI            = [%.1f, %.1f]  (robust SE %.2f)\n", w_lo, w_hi, se_wald))
cat(sprintf("Anderson-Rubin 95%% CI  = [%.1f, %.1f]  (grid inversion)\n", ar_lo, ar_hi))
cat(sprintf("AR set is a bounded interval: %s\n", ifelse(all(diff(match(inset,grid))==1),"yes","NO (unbounded/disjoint)")))
