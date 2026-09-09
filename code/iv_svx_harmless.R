## Harmless design (bias = 0, so H0* TRUE), nonlinear propensity swept over kappa.
## Question: does Directed-X hold 5% size as collinearity rises, or distort?
## Directed-S tests the SUFFICIENT condition (theta=0), which is FALSE here
## (curvature present) -> it should false-alarm regardless. This is the size/
## robustness face of "report both".
set.seed(20260902)
gen <- function(n, kappa, pi0=0.2, dp=0.4){
  X<-runif(n,-2,2); g<-X^2-4/3; e<-plogis(kappa*g); Z<-rbinom(n,1,e)
  pi1<-pi0+dp; V<-runif(n); D<-ifelse(Z==1,as.integer(V<=pi1),as.integer(V<=pi0))
  tau<-1; mu<-0                       # HARMLESS: const effect + linear (zero) level
  Y<-mu+tau*D+rnorm(n); list(X=X,Z=Z,D=D,Y=Y)
}
run<-function(d,degree=3,B=99){
  X<-d$X;Z<-d$Z;D<-d$D;Y<-d$Y;n<-length(Z);P<-poly(X,degree);M<-cbind(1,P);qi<-2:degree
  fit<-function(idx){Mi<-M[idx,,drop=FALSE];b<-qr.solve(Mi,Z[idx]);eh<-as.numeric(Mi%*%b)
    hh<-as.numeric(P[idx,qi,drop=FALSE]%*%b[1+qi]);zc<-Z[idx]-eh
    list(h=hh,Lbar=sum(zc*Y[idx])/sum(zc*D[idx]))}
  f0<-fit(seq_len(n));thD<-mean(f0$h*D);thY<-mean(f0$h*Y);gX<-mean(f0$h*(Y-f0$Lbar*D))
  bt<-matrix(NA,B,3)
  for(k in seq_len(B)){idx<-sample.int(n,n,TRUE);fk<-fit(idx)
    bt[k,]<-c(mean(fk$h*D[idx]),mean(fk$h*Y[idx]),mean(fk$h*(Y[idx]-fk$Lbar*D[idx])))}
  th<-c(thD,thY);W_S<-as.numeric(t(th)%*%solve(cov(bt[,1:2]),th))
  W_X<-gX^2/var(bt[,3]);c(DirS=1-pchisq(W_S,2),DirX=1-pchisq(W_X,1))
}
own<-function(kappa,N=200000){X<-runif(N,-2,2);g<-X^2-4/3;e<-plogis(kappa*g);mean(e*(1-e))}
kappas<-c(0.5,1,2,4,8,16)
tab<-t(sapply(seq_along(kappas),function(i){set.seed(200+i)
  MM<-t(replicate(300,run(gen(3000,kappas[i]))))
  c(kappa=kappas[i],ownvar_ZgivenX=own(kappas[i]),
    Directed_S_size=mean(MM[,"DirS"]<.05),Directed_X_size=mean(MM[,"DirX"]<.05))}))
cat("\n=== HARMLESS design (H0* true, bias=0): rejection = size, 5% nominal ===\n")
print(round(as.data.frame(tab),3),row.names=FALSE)
