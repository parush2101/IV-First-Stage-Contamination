set.seed(20260902)
d <- tryCatch(foreign::read.dta("final5.dta"), error=function(e) haven::read_dta("final5.dta"))
## Angrist-Lavy analysis sample: drop tiny/huge, valid scores
d <- subset(d, !is.na(avgmath) & !is.na(classize) & !is.na(cohsize) & !is.na(tipuach) &
              classize>0 & cohsize>=5 & cohsize<=200)
maim <- d$cohsize/(floor((d$cohsize-1)/40)+1)           # Maimonides predicted class size
Z <- maim; D <- d$classize; Y <- d$avgmath
cat("n=",length(Z)," corr(Z,cohsize)=",round(cor(Z,d$cohsize),3),"\n")
robust_t <- function(fit,nm){X<-model.matrix(fit);u<-resid(fit);nn<-nrow(X);k<-ncol(X)
  br<-solve(crossprod(X));V<-br%*%(t(X)%*%(X*u^2))%*%br*(nn/(nn-k));coef(fit)[nm]/sqrt(V[nm,nm])}
diag_iv <- function(Y,D,Z,Xl,Bt,B=999){Xl<-as.matrix(Xl);Bt<-as.matrix(Bt);n<-length(Z)
  rz<-resid(lm(Z~Xl));rd<-resid(lm(D~Xl));ry<-resid(lm(Y~Xl));bn<-sum(rz*ry)/sum(rz*rd)
  Fc<-robust_t(lm(D~Z+Xl),"Z")^2
  core<-function(ii){fl<-lm(Z[ii]~Xl[ii,,drop=F]);ff<-lm(Z[ii]~Xl[ii,,drop=F]+Bt[ii,,drop=F])
    eh<-fitted(ff);hh<-fitted(ff)-fitted(fl);rzi<-resid(fl)
    C<-mean(rzi*D[ii]);thD<-mean(hh*D[ii]);thY<-mean(hh*Y[ii]);zc<-Z[ii]-eh
    Lb<-sum(zc*Y[ii])/sum(zc*D[ii]);c(thD=thD,thY=thY,share=thD/C,Lbar=Lb,gX=mean(hh*(Y[ii]-Lb*D[ii])))}
  pt<-core(1:n);bt<-t(replicate(B,core(sample.int(n,n,TRUE))))
  WS<-as.numeric(t(pt[c("thD","thY")])%*%solve(cov(bt[,c("thD","thY")]),pt[c("thD","thY")]))
  WX<-pt["gX"]^2/var(bt[,"gX"])
  cat(sprintf("  beta=%.4f  F=%.1f  rho=%.3f[%.2f,%.2f]  p_S=%.4f  p_X=%.4f  corr=%.4f[%.3f,%.3f]\n",
    bn,Fc,pt["share"],quantile(bt[,"share"],.025),quantile(bt[,"share"],.975),
    1-pchisq(WS,2),1-pchisq(WX,1),pt["Lbar"],quantile(bt[,"Lbar"],.025),quantile(bt[,"Lbar"],.975)))}
cat("Maimonides, math score. Control: %disadvantaged + enrollment LINEAR; curvature: enroll^2, enroll^3\n")
Xl<-cbind(d$tipuach,d$cohsize); Bt<-resid(lm(cbind(d$cohsize^2,d$cohsize^3)~Xl))
diag_iv(Y,D,Z,Xl,Bt)
cat("Same, verbal score:\n")
diag_iv(d$avgverb,D,Z,Xl,Bt)
