set.seed(20260902)
d<-foreign::read.dta("401ksubs.dta")
robust_t <- function(fit,nm){X<-model.matrix(fit);u<-resid(fit);nn<-nrow(X);k<-ncol(X)
  br<-solve(crossprod(X));V<-br%*%(t(X)%*%(X*u^2))%*%br*(nn/(nn-k));coef(fit)[nm]/sqrt(V[nm,nm])}
diag_iv <- function(nm,Y,D,Z,Xl,Bt,B=1499){Xl<-as.matrix(Xl);Bt<-as.matrix(Bt);n<-length(Z)
  rz<-resid(lm(Z~Xl));rd<-resid(lm(D~Xl));ry<-resid(lm(Y~Xl));bn<-sum(rz*ry)/sum(rz*rd)
  Fc<-robust_t(lm(D~Z+Xl),"Z")^2
  core<-function(ii){fl<-lm(Z[ii]~Xl[ii,,drop=F]);ff<-lm(Z[ii]~Xl[ii,,drop=F]+Bt[ii,,drop=F])
    eh<-fitted(ff);hh<-fitted(ff)-fitted(fl);rzi<-resid(fl)
    C<-mean(rzi*D[ii]);zc<-Z[ii]-eh;Lb<-sum(zc*Y[ii])/sum(zc*D[ii])
    c(thD=mean(hh*D[ii]),thY=mean(hh*Y[ii]),share=mean(hh*D[ii])/C,Lbar=Lb,gX=mean(hh*(Y[ii]-Lb*D[ii])))}
  pt<-core(1:n);bt<-t(replicate(B,core(sample.int(n,n,TRUE))))
  WS<-as.numeric(t(pt[c("thD","thY")])%*%solve(cov(bt[,c("thD","thY")]),pt[c("thD","thY")]));WX<-pt["gX"]^2/var(bt[,"gX"])
  flag<-ifelse(1-pchisq(WX,1)<0.05,"  <<< BIAS FLAG","")
  cat(sprintf("%-26s n=%d beta=%8.2f F=%7.0f rho=%5.2f p_S=%.4f p_X=%.4f corr=%8.2f[%.1f,%.1f]%s\n",
    nm,n,bn,Fc,pt["share"],1-pchisq(WS,2),1-pchisq(WX,1),pt["Lbar"],quantile(bt[,"Lbar"],.025),quantile(bt[,"Lbar"],.975),flag))
  invisible(pt)}
# reference: flexible-income 2SLS (income splines) = the "correct" answer
library(splines)
fx<-lm(nettfa~p401k+ns(inc,5)+ns(age,4)+marr+male+fsize,d)  # not IV, just for context
cat("--- 401(k): effect of participation (p401k) on net financial assets, instrument e401k ---\n")
cat("A) income & age controlled LINEARLY (common default):\n  ")
Xl<-with(d,cbind(inc,age,marr,male,fsize)); Bt<-resid(lm(with(d,cbind(inc^2,inc^3,age^2,inc*age))~Xl))
diag_iv("401k linear income",d$nettfa,d$p401k,d$e401k,Xl,Bt)
cat("B) income quadratic (incsq), age quadratic:\n  ")
Xl<-with(d,cbind(inc,incsq,age,agesq,marr,male,fsize)); Bt<-resid(lm(with(d,cbind(inc^3,inc^4,age^3,inc*age))~Xl))
diag_iv("401k quadratic income",d$nettfa,d$p401k,d$e401k,Xl,Bt)
cat("C) robustness of the flag: different curvature basis (splines of income) under linear control:\n  ")
Xl<-with(d,cbind(inc,age,marr,male,fsize)); Bt<-resid(lm(ns(d$inc,6)~Xl))
diag_iv("401k linear, spline basis",d$nettfa,d$p401k,d$e401k,Xl,Bt)
