## Final honest empirical panel for the paper. Reports ALL designs.
## (CollegeDistance excluded: its 'wage' is a 41-value regional average, not an
##  individual outcome, so it is not a valid returns-to-schooling IV.)
suppressMessages(library(AER)); suppressMessages(library(hdm)); set.seed(20260902)
robust_t <- function(fit,nm){X<-model.matrix(fit);u<-resid(fit);nn<-nrow(X);k<-ncol(X)
  br<-solve(crossprod(X));V<-br%*%(t(X)%*%(X*u^2))%*%br*(nn/(nn-k));coef(fit)[nm]/sqrt(V[nm,nm])}
diag_iv <- function(nm,Y,D,Z,Xl,Bt,B=999){Xl<-as.matrix(Xl);Bt<-as.matrix(Bt);n<-length(Z)
  rz<-resid(lm(Z~Xl));rd<-resid(lm(D~Xl));ry<-resid(lm(Y~Xl));bn<-sum(rz*ry)/sum(rz*rd)
  Fc<-robust_t(lm(D~Z+Xl),"Z")^2
  core<-function(ii){fl<-lm(Z[ii]~Xl[ii,,drop=F]);ff<-lm(Z[ii]~Xl[ii,,drop=F]+Bt[ii,,drop=F])
    eh<-fitted(ff);hh<-fitted(ff)-fitted(fl);rzi<-resid(fl)
    C<-mean(rzi*D[ii]);zc<-Z[ii]-eh;Lb<-sum(zc*Y[ii])/sum(zc*D[ii])
    c(thD=mean(hh*D[ii]),thY=mean(hh*Y[ii]),share=mean(hh*D[ii])/C,Lbar=Lb,gX=mean(hh*(Y[ii]-Lb*D[ii])))}
  pt<-core(1:n);bt<-t(replicate(B,core(sample.int(n,n,TRUE))))
  WS<-as.numeric(t(pt[c("thD","thY")])%*%solve(cov(bt[,c("thD","thY")]),pt[c("thD","thY")]));WX<-pt["gX"]^2/var(bt[,"gX"])
  cat(sprintf("%-18s n=%4d beta=%8.4f F=%6.1f rho=%5.2f p_S=%.3f p_X=%.3f corr=%8.4f\n",
    nm,n,bn,Fc,pt["share"],1-pchisq(WS,2),1-pchisq(WX,1),pt["Lbar"]))}
# Card
d<-foreign::read.dta("card.dta");ct<-c("exper","expersq","black","south","smsa","smsa66",paste0("reg66",2:9))
d<-d[complete.cases(d[,c("lwage","educ","nearc4",ct)]),];Xl<-model.matrix(~.,d[,ct])[,-1]
Bt<-resid(lm(with(d,cbind(exper^3,exper*black,exper*south,exper*smsa,black*south,black*smsa,south*smsa,exper*smsa66))~Xl))
diag_iv("Card (schooling)",d$lwage,d$educ,d$nearc4,Xl,Bt)
# Mroz
data("PSID1976");m<-PSID1976[PSID1976$participation=="yes",];Xl<-model.matrix(~experience+I(experience^2)+age+youngkids+oldkids+city,m)[,-1]
Bt<-resid(lm(with(m,cbind(experience^3,age^2,age^3,experience*age,experience*youngkids,age*(city=="yes")))~Xl))
diag_iv("Mroz (schooling)",log(m$wage),m$education,m$meducation,Xl,Bt)
# CigarettesSW
data("CigarettesSW");cs<-subset(CigarettesSW,year=="1995");ri<-log(cs$income/cs$population/cs$cpi)
Xl<-cbind(ri);Bt<-resid(lm(cbind(ri^2,ri^3)~Xl));diag_iv("Cigarettes (demand)",log(cs$packs),log(cs$price/cs$cpi),cs$taxs/cs$cpi,Xl,Bt)
# AJR
data("AJR");Xl<-cbind(AJR$Latitude);Bt<-resid(lm(cbind(AJR$Latitude^2,AJR$Latitude^3)~Xl))
diag_iv("AJR (institutions)",AJR$GDP,AJR$Exprop,AJR$logMort,Xl,Bt)
# Maimonides
d<-tryCatch(foreign::read.dta("final5.dta"),error=function(e)haven::read_dta("final5.dta"))
d<-subset(d,!is.na(avgmath)&!is.na(classize)&!is.na(cohsize)&!is.na(tipuach)&classize>0&cohsize>=5&cohsize<=200)
Z<-d$cohsize/(floor((d$cohsize-1)/40)+1);Xl<-cbind(d$tipuach,d$cohsize);Bt<-resid(lm(cbind(d$cohsize^2,d$cohsize^3)~Xl))
diag_iv("Maimonides (math)",d$avgmath,d$classize,Z,Xl,Bt)
diag_iv("Maimonides (verbal)",d$avgverb,d$classize,Z,Xl,Bt)
