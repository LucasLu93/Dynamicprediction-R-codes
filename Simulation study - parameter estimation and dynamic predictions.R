library(INLA)
library(mvtnorm)
library(countreg)
library(IDPmisc)
library(distributions3)
library(dplyr)
inla.setOption(num.threads=6) 
set.seed(1)

######### True values for parameters ###########################################
# random effects
sd.Y.R<-1.5 # related to recurrent events
var.Y.R<-sd.Y.R^2
sd.Y.L<-1.8 # related to longitudinal outcomes
var.Y.L<-sd.Y.L^2
cor.Y<-0.9
cov.Y<-cor.Y*sqrt(var.Y.R*var.Y.L)
vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),2,2)

sd.Z.R<-0.5 # related to recurrent events
var.Z.R<-sd.Z.R^2
sd.Z.L<-0.7 # related to longitudinal outcomes
var.Z.L<-sd.Z.L^2
cor.Z<-0.3
cov.Z<-cor.Z*sqrt(var.Z.R*var.Z.L)
vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),2,2)

# terminal event
beta<-c(-5.5,-1.1)
xi<--0.5
lambda.d<-0.7 # shape for weibull distribution 
rho.d<-33.5 # scale for weibull distribution

# recurrent event
gamma<-c(0.3,-0.3)
alpha<-0.4
lambda.r<-1.2 # shape for weibull distribution
rho.r<-4.6 # scale for weibull distribution 

# longitudinal outcomes
chi<-c(1.6,-0.3,-0.2)
theta<-c(-4.9,1.6,-1.1)
omega.l<--0.6
omega.z<-0.7
tau.l<-1.2 
tau.z<-1.2 
varphi<-25.5 # dispersion parameter 

######### Generating simulated data ############################################
n2<-200 # number of subjects
n3<-20 # number of families
m<-10 # maximum number of observed times for each subjects    
id<-1:n2
id.2<-rep(id,each=m)

famid<-rep(1:n3,each=10)
famid.2<-rep(famid,each=m)

# random effects
Y<-rmvnorm(n2,c(0,0),vcov.Y)
Z<-rmvnorm(n3,c(0,0),vcov.Z)

Y.2<-Y[id.2,] # subject-specific random effect
Z.2<-Z[famid.2,] # family-specific random effect

# terminal event
x1<-rbinom(n2,1,0.4) 
x1.d<-x1
x2<-rnorm(n3,-0.2,0.5) 
x2.d<-x2[famid] 

mu.d<-apply(cbind(x1.d,x2.d,Y[,2],Z[famid,][,2],Y[,1],Z[famid,][,1]),1,function(x){
  beta[1]*x[1]+beta[2]*x[2]+x[3]+x[4]+x[5]+xi*x[6]})
surv.d<-sapply(mu.d,function(x){rweibull(1,lambda.d,rho.d/exp(x))}) 

q<-rexp(n2,0.1)
status.d<-rep(1,n2) 
for (i in 1:n2){
  if (surv.d[i]>q[i]) {
    surv.d[i]<-q[i]
    status.d[i]<-0
  }
}

x1.d<-rep(x1.d,each=m)
x2.d<-rep(x2.d,each=m)
surv.d<-rep(surv.d,each=m)
status.d<-rep(status.d,each=m)

# recurrent event
x1.r<-x1[id.2]
x2.r<-x2[famid.2] 

mu.r<-apply(cbind(x1.r,x2.r,Y.2[,1],Z.2[,1]),1,function(x){
  gamma[1]*x[1]+gamma[2]*x[2]+alpha*x[3]+x[4]})
surv.r<-sapply(mu.r,function(x){rweibull(1,lambda.r,rho.r/exp(x))})

cum.gap<-surv.r[1] 
for (i in 1:n2){
  for (j in 1:m){
    if (j==1) (cum.gap[(i-1)*m+j]<-surv.r[(i-1)*m+j])
    if (j>1) (cum.gap[(i-1)*m+j]<-cum.gap[(i-1)*m+j-1]+surv.r[(i-1)*m+j])
  }
}

status.r<-rep(1,n2*m) 
v<-rep(1,n2*m) 
for (i in 1:n2){
  for (j in 1:m){
    if (j==1 & cum.gap[(i-1)*m+1]>surv.d[(i-1)*m+1]) {
      v[((i-1)*m+2):((i-1)*m+m)]<-0
      surv.r[(i-1)*m+1]<-surv.d[(i-1)*m+1]
      status.r[(i-1)*m+1]<-0
    }
    if (j>1 & j<m & cum.gap[(i-1)*m+j]>surv.d[(i-1)*m+j]) {
      v[((i-1)*m+j+1):((i-1)*m+m)]<-0
      surv.r[(i-1)*m+j]<-surv.d[(i-1)*m+j]-cum.gap[(i-1)*m+j-1]
      status.r[(i-1)*m+j]<-0
    }
    if (j==m & cum.gap[(i-1)*m+m]>surv.d[(i-1)*m+m]) {
      surv.r[(i-1)*m+m]<-surv.d[(i-1)*m+m]-cum.gap[(i-1)*m+m-1]
      status.r[(i-1)*m+m]<-0
    }
  }
}

# longitudinal outcomes
T.R<-surv.r 
log.T.R<-log(T.R) # offset term
log.T.R[is.na(log.T.R)]<-0
x1.l<-x1[id.2]
x2.l<-x2[famid.2] 

# logistic component
mu.l<-apply(cbind(log.T.R,x1.l,x2.l,Y.2[,2],Z.2[,2]),1,function(x){
  x[1]+chi[1]+chi[2]*x[2]+chi[3]*x[3]+omega.l*x[4]+tau.l*x[5]}) 
z<-sapply(mu.l,function(x){rbinom(1,1,exp(x)/(1+exp(x)))})

# negative binomial component
mu.z<-apply(cbind(log.T.R,x1.l,x2.l,Y.2[,2],Z.2[,2]),1,function(x){
  x[1]+theta[1]+theta[2]*x[2]+theta[3]*x[3]+omega.z*x[4]+tau.z*x[5]}) 

y<-NULL 
for (i in 1:(n2*m)){
  if (z[i]==1) (y[i]<-0)
  if (z[i]==0) (y[i]<-rztnbinom(1,mu=exp(mu.z[i]),theta=varphi))
}

# generating simulated data
data<-cbind(id.2,famid.2,T.R,x1.l,x2.l,x1.r,x2.r,surv.r,status.r,x1.d,x2.d,surv.d,status.d,y)
data<-data[v==1,]
d<-NULL
n<-nrow(data)
for (i in 1:(n-1)){
  if (data[i,1]==data[i+1,1]) {
    d[i]<-0} else {
      d[i]<-1}}
d[n]<-1
data<-data.frame(cbind(data,d))
data$cumgap<-ave(data$T.R,data$id.2,FUN=cumsum)
data$n<-ave(data$id.2,data$id.2,FUN=length)
data$time<-ave(data$id.2,data$id.2,FUN=seq_along)+1
data$visit<-1-data$d

######### Fitting proposed trivariate joint model ##############################
# preparation for longitudinal outcomes
data.l<-data[data$y!=0,]
n4<-nrow(data.l)

# preparation for recurrent events
time.R<-data$surv.r 
delta.R.fij<-data$status.r 

# preparation for terminal event
data.d<-data[data$d==1,]
time.D<-data.d$surv.d 
delta.D.fi<-data.d$status.d 

# other preparations
n1<-nrow(data)
z<-rep(0,n1)
z[data$y==0]<-1

l.long <- c(z, rep(NA, n4),  rep(NA, n1), rep(NA, n2))
z.long <- c(rep(NA, n1), data.l$y, rep(NA, n1), rep(NA, n2))
y.recu <- inla.surv(time = c(rep(NA, n1+n4), time.R, rep(NA, n2)), event = c(rep(NA,n1+n4), delta.R.fij, rep(NA, n2)))
y.term <- inla.surv(time = c(rep(NA, n1+n4+n1), time.D), event = c(rep(NA, n1+n4+n1), delta.D.fi))
y.joint <- list(l.long, z.long, y.recu, y.term)

linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n4),rep(1,n1),rep(2,n2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                               l.T.R = c(log(data$T.R), rep(0, n1+n4), rep(0, n2)),  # related to logit(pi): logistic component for longitudinal outcomes
                               intercept.l = c(rep(1,n1), rep(0, n1+n4), rep(0,n2)),
                               x1.l = c(data$x1.l, rep(0, n1+n4), rep(0, n2)), 
                               x2.l = c(data$x2.l, rep(0, n1+n4), rep(0, n2)),
                               z.T.R = c(rep(0, n1), log(data.l$T.R), rep(0,n1), rep(0, n2)),  # related to log(mu): negative binomial component for longitudinal outcomes
                               intercept.z = c(rep(0,n1), rep(1,n4), rep(0,n1), rep(0,n2)),
                               x1.z = c(rep(0, n1), data.l$x1.l, rep(0,n1), rep(0, n2)),
                               x2.z = c(rep(0, n1), data.l$x2.l, rep(0,n1), rep(0, n2)),
                               x1.r = c(rep(0, n1+n4), data$x1.r, rep(0, n2)), # relative to hazard function for recurrent events
                               x2.r = c(rep(0, n1+n4), data$x2.r, rep(0, n2)),
                               x1.d = c(rep(0, n1+n4+n1), data.d$x1.d), # related to hazard function for terminal event
                               x2.d = c(rep(0, n1+n4+n1), data.d$x2.d))

random.covariate <- list(l.YL = c(data$id+n2, rep(NA, n1+n4), rep(NA, n2)), 
                         l.ZL = c(data$famid+n3, rep(NA, n1+n4), rep(NA, n2)), 
                         z.YL = c(rep(NA, n1), data.l$id+n2, rep(NA, n1), rep(NA, n2)), 
                         z.ZL = c(rep(NA, n1), data.l$famid+n3, rep(NA, n1), rep(NA, n2)), 
                         r.YR = c(rep(NA, n1+n4), data$id, rep(NA, n2)), 
                         r.ZR = c(rep(NA, n1+n4), data$famid, rep(NA, n2)), 
                         d.YR = c(rep(NA, n1+n4+n1), 1:n2), 
                         d.ZR = c(rep(NA, n1+n4+n1), data.d$famid),
                         d.YL = c(rep(NA, n1+n4+n1), (1:n2)+n2), 
                         d.ZL = c(rep(NA, n1+n4+n1), data.d$famid+n3))

data.f <- c(linear.covariate,random.covariate)
data.f$Y <- y.joint

formula = Y ~ - 1 + mu + intercept.l + x1.l + x2.l + offset(l.T.R) +
  intercept.z + x1.z + x2.z + offset(z.T.R) +
  x1.r + x2.r + 
  x1.d + x2.d + 
  f(d.YR, model="iidkd", order=2, n=2*n2, hyper = list(theta1 = list(param = c(5, 1, 1, 0)))) +
  f(d.YL, copy="d.YR") +
  f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(l.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(r.ZR, model="iidkd", order=2, n=2*n3, hyper = list(theta1 = list(param = c(5, 1, 1, 0)))) +
  f(d.ZL, copy="r.ZR") +
  f(l.ZL, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.ZL, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(d.ZR, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

inla.model<-inla(formula, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),data = data.f,
                 control.compute=list(dic=FALSE,cpo=FALSE,waic=FALSE),
                 control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                 control.inla = list(int.strategy = "eb"), inla.mode = "compact", safe = TRUE,
                 control.fixed = list(prec.intercept = 0.1))

inla.model.result<-summary(inla.model)

# parameter estimates
# fixed effects
est.fixed<-inla.model.result$fixed[,1][-c(1,2)]

# hyperparameters
# scale parameter
r.r.mean<-inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model$marginals.fixed[[1]]), silent = TRUE)$mean
r.d.mean<-inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model$marginals.fixed[[2]]), silent = TRUE)$mean

# dispersion
d.mean<-inla.zmarginal(inla.tmarginal(function(x) exp(x),inla.model$internal.marginals.hyperpar[[1]]), silent = TRUE)$mean

# standard deviation and correlation coefficients of subject-specific random effects Y
mcsamples <- inla.iidkd.sample(10^4, inla.model, "d.YR", return.cov=FALSE)
sdcor.Y <- matrix(unlist(mcsamples), nrow = 2^2)
sdcor.Y.mean <- rowMeans(sdcor.Y)[-3]

# standard deviation and correlation coefficients of family-specific random effects Z
mcsamples <- inla.iidkd.sample(10^4, inla.model, "r.ZR", return.cov=FALSE)
sdcor.Z <- matrix(unlist(mcsamples), nrow = 2^2)
sdcor.Z.mean <- rowMeans(sdcor.Z)[-3]

est.hyper<-c(r.r.mean,r.d.mean,d.mean,inla.model.result$hyperpar[,1][c(2,3)],sdcor.Y.mean,sdcor.Z.mean,inla.model.result$hyperpar[,1][seq(10,15)])

# generate results
est.total<-c(est.fixed,est.hyper)
names(est.total)<-c("intercept.l","x1.l","x2.l","intercept.z","x1.z","x2.z","x1.r","x2.r","x1.d","x2.d",
                       "rho.r","rho.d","varphi","lambda.r","lambda.d","sd.Y.R","cor.Y","sd.Y.L","sd.Z.R","cor.Z","sd.Z.L",
                       "alpha","omega.l","omega.z","tau.l","tau.z","xi")

######### Generating dynamic predictions for terminal event ####################
eta.te12<-function(para1,i1,t1,data){
  param<-as.numeric(para1)
  
  # longitudinal outcomes
  # logistic component
  chi0<-param[1]
  chi1<-param[2]
  chi2<-param[3]
  omega.l<-param[23]
  tau.l<-param[25]
  dispersion<-param[13]
  
  # negative binomial component
  theta0<-param[4]
  theta1<-param[5]
  theta2<-param[6]
  omega.z<-param[24]
  tau.z<-param[26]
  
  # recurrent events
  gamma1<-param[7]
  gamma2<-param[8]
  a.R<-param[14]
  b.R<-param[11]
  alpha<-param[22]
  
  # terminal event
  beta1<-param[9]
  beta2<-param[10]
  a.D<-param[15]
  b.D<-param[12]
  xi<-param[27]
  
  # random effects
  var.Y.R<-param[16]^2
  var.Y.L<-param[18]^2
  rho.Y<-param[17]
  cov.Y<-rho.Y*sqrt(var.Y.R)*sqrt(var.Y.L)
  vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),ncol=2)
  W.Y<-solve(vcov.Y)
  R.Y<-chol(W.Y)
  L.Y<-t(R.Y)
  var.Z.R<-param[19]^2
  var.Z.L<-param[21]^2
  rho.Z<-param[20]
  cov.Z<-rho.Z*sqrt(var.Z.R)*sqrt(var.Z.L)
  vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),ncol=2)
  W.Z<-solve(vcov.Z)
  R.Z<-chol(W.Z)
  L.Z<-t(R.Z)
  
  # family members' information
  famid_select <- data$famid.2[data$id.2==i1][1]
  fam_hist <- data[data$famid.2==famid_select, ]
  fam_hist$delta.R <- ifelse(fam_hist$cumgap <= t1,1,0)
  jstar.i <- aggregate(fam_hist$delta.R, by=list(fam_hist$id.2), sum)[,2]+1
  n.i <- fam_hist$n[fam_hist$time==2]
  fam_hist$jstar <- rep(jstar.i, n.i)
  
  # preparation for terminal event 
  data.term<-fam_hist[fam_hist$d==1, ]
  X.D.f<-with(data.term,cbind(x1.d,x2.d))
  id.D.f<-data.term$id.2
  famid.D.f<-data.term$famid.2
  Delta.D.f<-data.term$status.d
  T.D.f<-data.term$surv.d
  T.fh.D <- pmin(T.D.f,t1)
  delta.fh.D <- ifelse(T.D.f <= t1, Delta.D.f, 0)
  
  # preparation for recurrent events
  data.rec1<-fam_hist[fam_hist$jstar==1 & fam_hist$time==2,]
  X.R.f1<-with(data.rec1,cbind(x1.r,x2.r))
  id.R.f1<-data.rec1$id.2
  famid.R.f1<-data.rec1$famid.2
  Delta.R.f1<-data.rec1$delta.R
  T.R.f1<-rep(t1,nrow(data.rec1))
  data.rec2<-fam_hist[fam_hist$jstar>1 & fam_hist$time <= fam_hist$jstar+1,] 
  X.R.f2 <- with(data.rec2,cbind(x1.r,x2.r))
  id.R.f2<-data.rec2$id.2
  famid.R.f2<-data.rec2$famid.2
  Delta.R.f2<-data.rec2$delta.R*data.rec2$visit
  T.R.f2<-with(data.rec2, ifelse(delta.R==1, T.R, t1 - cumgap + T.R))
  T.R.f2[T.R.f2<1e-10]<-0
  data.rec<-rbind(data.rec1,data.rec2)
  X.R.f<-rbind(X.R.f1,X.R.f2)
  id.R.f<-c(id.R.f1,id.R.f2)
  famid.R.f<-c(famid.R.f1,famid.R.f2)
  Delta.R.f<-c(Delta.R.f1,Delta.R.f2)
  T.R.f<-c(T.R.f1,T.R.f2)
  
  # preparation for longitudinal outcomes
  # logistic component
  data.long<-data.rec2[data.rec2$delta.R==1,]
  X.L.f<-with(data.long,cbind(x1.l,x2.l))
  z<-rep(0,nrow(data.long))
  z[data.long$y==0]<-1
  
  # negative binomial component
  data.long.z<-data.long[data.long$y != 0,]
  X.L.f.z<-with(data.long.z,cbind(x1.l,x2.l))
  
  # other preparations
  n1<-nrow(data.long)
  n2<-nrow(data.long.z)
  n3<-nrow(data.rec)
  n4<-nrow(data.term)
  
  # construct model
  l.long<-c(z,rep(NA,n2+n3+n4))
  z.long <- c(rep(NA,n1),data.long.z$y,rep(NA,n3+n4))
  y.recu<-INLA::inla.surv(time=c(rep(NA,n1+n2),T.R.f,rep(NA,n4)), event=c(rep(NA,n1+n2),Delta.R.f,rep(NA,n4)))
  y.term<-INLA::inla.surv(time=c(rep(NA,n1+n2+n3),T.fh.D), event=c(rep(NA,n1+n2+n3),delta.fh.D))
  y.joint <- list(l.long, z.long, y.recu, y.term)
  
  linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n2),rep(1,n3),rep(2,n4))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.Time = c(log(data.long$T.R), rep(NA,n2+n3+n4)), # related to logit(pi): logistic component for longitudinal outcomes 
                                 l.x1 = c(data.long$x1.l, rep(NA,n2+n3+n4)),
                                 l.x2 = c(data.long$x2.l, rep(NA,n2+n3+n4)),
                                 z.Time = c(rep(NA,n1), log(data.long.z$T.R), rep(NA,n3+n4)), # related to log(mu): negative binomial component for longitudinal outcomes
                                 z.x1 = c(rep(NA,n1), data.long.z$x1.l, rep(NA,n3+n4)),
                                 z.x2 = c(rep(NA,n1), data.long.z$x2.l, rep(NA,n3+n4)),
                                 r.x1 = c(rep(NA,n1+n2), data.rec$x1.r, rep(NA,n4)), # relative to hazard function for recurrent events
                                 r.x2 = c(rep(NA,n1+n2), data.rec$x2.r, rep(NA,n4)),
                                 d.x1 = c(rep(NA,n1+n2+n3), data.term$x1.d), # related to hazard function for terminal event
                                 d.x2 = c(rep(NA,n1+n2+n3), data.term$x2.d))
  
  random.covariate <- list(l.YL = c(data.long$id.2+200, rep(NA,n2+n3+n4)), 
                           l.ZL = c(data.long$famid.2+20, rep(NA,n2+n3+n4)), 
                           z.YL = c(rep(NA,n1), data.long.z$id.2+200, rep(NA,n3+n4)), 
                           z.ZL = c(rep(NA,n1), data.long.z$famid.2+20, rep(NA,n3+n4)), 
                           r.YR = c(rep(NA,n1+n2), data.rec$id.2, rep(NA,n4)), 
                           r.ZR = c(rep(NA,n1+n2), data.rec$famid.2, rep(NA,n4)), 
                           d.YR = c(rep(NA,n1+n2+n3), data.term$id.2), 
                           d.ZR = c(rep(NA,n1+n2+n3), data.term$famid.2), 
                           d.YL = c(rep(NA,n1+n2+n3), data.term$id.2+200), 
                           d.ZL = c(rep(NA,n1+n2+n3), data.term$famid.2+20))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  # Trivariate joint model
  formula1 = Y ~ -1 + offset(c(rep(chi0,n1),rep(theta0,n2),rep(-log(b.R),n3),rep(-log(b.D),n4))) +
    l.x1 + l.x2 + offset(l.Time) +
    z.x1 + z.x2 + offset(z.Time) +
    r.x1 + r.x2 +
    d.x1 + d.x2 +
    f(d.YR, model="iidkd", order=2, n=2*200, hyper = list(theta1 = list(initial=log(L.Y[1,1]),fixed=TRUE),
                                                         theta2 = list(initial=log(L.Y[2,2]),fixed=TRUE),
                                                         theta3 = list(initial=L.Y[2,1],fixed=TRUE))) +
    f(d.YL, copy="d.YR") +
    f(r.YR, copy="d.YR", hyper = list(beta = list(initial=alpha, fixed=TRUE))) +
    f(l.YL, copy="d.YR", hyper = list(beta = list(initial=omega.l, fixed=TRUE)), n = length(sort(unique(l.YL)))) +
    f(z.YL, copy="d.YR", hyper = list(beta = list(initial=omega.z, fixed=TRUE)), n = length(sort(unique(z.YL)))) +
    f(r.ZR, model="iidkd", order=2, n=2*20, hyper = list(theta1 = list(initial=log(L.Z[1,1]),fixed=TRUE),
                                                        theta2 = list(initial=log(L.Z[2,2]),fixed=TRUE),
                                                        theta3 = list(initial=L.Z[2,1],fixed=TRUE))) +
    f(d.ZL, copy="r.ZR") +
    f(l.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.l, fixed=TRUE)), n = length(sort(unique(l.ZL)))) +
    f(z.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.z, fixed=TRUE)), n = length(sort(unique(z.ZL)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(initial=xi, fixed=TRUE))) 
  
  model1 <- INLA::inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                       data = data.f, control.compute=list(dic=FALSE,cpo=FALSE,waic=FALSE, return.marginals.predictor = FALSE,config = FALSE),
                       control.family = list(list(),
                                             list(hyper=list(size=list(initial=log(dispersion),fixed=TRUE),
                                                             prob=list(initial=-10,fixed = TRUE))),
                                             list(variant=1,hyper=list(alpha=list(initial=log(a.R),fixed=TRUE))),
                                             list(variant=1,hyper=list(alpha=list(initial=log(a.D),fixed=TRUE)))),
                       control.inla = list(int.strategy = "eb"), inla.mode = "compact", safe = TRUE, 
                       control.fixed = list(mean=list(l.x1=chi1,l.x2=chi2,
                                                      z.x1=theta1,z.x2=theta2,
                                                      r.x1=gamma1,r.x2=gamma2,
                                                      d.x1=beta1,d.x2=beta2),
                                            prec=list(l.x1=1e10,l.x2=1e10,
                                                      z.x1=1e10,z.x2=1e10,
                                                      r.x1=1e10,r.x2=1e10,
                                                      d.x1=1e10,d.x2=1e10)),
                       control.predictor = list(link=1))
  
  # model 1
  eta.c1<-model1$summary.linear.predictor$mean[n1+n2+n3+c(1:n4)]
  
  # model 2
  eta.c2<-model1$summary.linear.predictor$mean[n1+n2+n3+c(1:n4)]-
    xi*model1$summary.random$d.ZL$mean[famid_select]-model1$summary.random$d.ZL$mean[famid_select+20]
  
  return(rbind(eta.c1,eta.c2))
}

eta.te3<-function(para1,i1,t1,data){
  param<-as.numeric(para1)
  
  # longitudinal outcomes
  # logistic component
  chi0<-param[1]
  chi1<-param[2]
  chi2<-param[3]
  omega.l<-param[23]
  tau.l<-param[25]
  dispersion<-param[13]
  
  # negative binomial component
  theta0<-param[4]
  theta1<-param[5]
  theta2<-param[6]
  omega.z<-param[24]
  tau.z<-param[26]
  
  # recurrent events
  gamma1<-param[7]
  gamma2<-param[8]
  a.R<-param[14]
  b.R<-param[11]
  alpha<-param[22]
  
  # terminal event
  beta1<-param[9]
  beta2<-param[10]
  a.D<-param[15]
  b.D<-param[12]
  xi<-param[27]
  
  # random effects
  var.Y.R<-param[16]^2
  var.Y.L<-param[18]^2
  rho.Y<-param[17]
  cov.Y<-rho.Y*sqrt(var.Y.R)*sqrt(var.Y.L)
  vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),ncol=2)
  W.Y<-solve(vcov.Y)
  R.Y<-chol(W.Y)
  L.Y<-t(R.Y)
  var.Z.R<-param[19]^2
  var.Z.L<-param[21]^2
  rho.Z<-param[20]
  cov.Z<-rho.Z*sqrt(var.Z.R)*sqrt(var.Z.L)
  vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),ncol=2)
  W.Z<-solve(vcov.Z)
  R.Z<-chol(W.Z)
  L.Z<-t(R.Z)
  
  # family members' information
  famid_select <- data$famid.2[data$id.2==i1][1]
  fam_hist <- data[data$famid.2==famid_select & data$id.2==i1, ]
  fam_hist$delta.R <- ifelse(fam_hist$cumgap <= t1,1,0)
  jstar.i <- aggregate(fam_hist$delta.R, by=list(fam_hist$id.2), sum)[,2]+1
  n.i <- fam_hist$n[fam_hist$time==2]
  fam_hist$jstar <- rep(jstar.i, n.i)
  
  # preparation for terminal event 
  data.term<-fam_hist[fam_hist$d==1, ]
  X.D.f<-with(data.term,cbind(x1.d,x2.d))
  id.D.f<-data.term$id.2
  famid.D.f<-data.term$famid.2
  Delta.D.f<-data.term$status.d
  T.D.f<-data.term$surv.d
  T.fh.D <- pmin(T.D.f,t1)
  delta.fh.D <- ifelse(T.D.f <= t1, Delta.D.f, 0)
  
  # preparation for recurrent events
  data.rec1<-fam_hist[fam_hist$jstar==1 & fam_hist$time==2,]
  X.R.f1<-with(data.rec1,cbind(x1.r,x2.r))
  id.R.f1<-data.rec1$id.2
  famid.R.f1<-data.rec1$famid.2
  Delta.R.f1<-data.rec1$delta.R
  T.R.f1<-rep(t1,nrow(data.rec1))
  data.rec2<-fam_hist[fam_hist$jstar>1 & fam_hist$time <= fam_hist$jstar+1,] 
  X.R.f2 <- with(data.rec2,cbind(x1.r,x2.r))
  id.R.f2<-data.rec2$id.2
  famid.R.f2<-data.rec2$famid.2
  Delta.R.f2<-data.rec2$delta.R*data.rec2$visit
  T.R.f2<-with(data.rec2, ifelse(delta.R==1, T.R, t1 - cumgap + T.R))
  T.R.f2[T.R.f2<1e-10]<-0
  data.rec<-rbind(data.rec1,data.rec2)
  X.R.f<-rbind(X.R.f1,X.R.f2)
  id.R.f<-c(id.R.f1,id.R.f2)
  famid.R.f<-c(famid.R.f1,famid.R.f2)
  Delta.R.f<-c(Delta.R.f1,Delta.R.f2)
  T.R.f<-c(T.R.f1,T.R.f2)
  
  # preparation for longitudinal outcomes
  # logistic component
  data.long<-data.rec2[data.rec2$delta.R==1,]
  X.L.f<-with(data.long,cbind(x1.l,x2.l))
  z<-rep(0,nrow(data.long))
  z[data.long$y==0]<-1
  
  # negative binomial component
  data.long.z<-data.long[data.long$y != 0,]
  X.L.f.z<-with(data.long.z,cbind(x1.l,x2.l))
  
  # other preparations
  n1<-nrow(data.long)
  n2<-nrow(data.long.z)
  n3<-nrow(data.rec)
  n4<-nrow(data.term)
  
  # construct model
  l.long<-c(z,rep(NA,n2+n3+n4))
  z.long <- c(rep(NA,n1),data.long.z$y,rep(NA,n3+n4))
  y.recu<-INLA::inla.surv(time=c(rep(NA,n1+n2),T.R.f,rep(NA,n4)), event=c(rep(NA,n1+n2),Delta.R.f,rep(NA,n4)))
  y.term<-INLA::inla.surv(time=c(rep(NA,n1+n2+n3),T.fh.D), event=c(rep(NA,n1+n2+n3),delta.fh.D))
  y.joint <- list(l.long, z.long, y.recu, y.term)
  
  linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n2),rep(1,n3),rep(2,n4))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.Time = c(log(data.long$T.R), rep(NA,n2+n3+n4)), # related to logit(pi): logistic component for longitudinal outcomes 
                                 l.x1 = c(data.long$x1.l, rep(NA,n2+n3+n4)),
                                 l.x2 = c(data.long$x2.l, rep(NA,n2+n3+n4)),
                                 z.Time = c(rep(NA,n1), log(data.long.z$T.R), rep(NA,n3+n4)), # related to log(mu): negative binomial component for longitudinal outcomes
                                 z.x1 = c(rep(NA,n1), data.long.z$x1.l, rep(NA,n3+n4)),
                                 z.x2 = c(rep(NA,n1), data.long.z$x2.l, rep(NA,n3+n4)),
                                 r.x1 = c(rep(NA,n1+n2), data.rec$x1.r, rep(NA,n4)), # relative to hazard function for recurrent events
                                 r.x2 = c(rep(NA,n1+n2), data.rec$x2.r, rep(NA,n4)),
                                 d.x1 = c(rep(NA,n1+n2+n3), data.term$x1.d), # related to hazard function for terminal event
                                 d.x2 = c(rep(NA,n1+n2+n3), data.term$x2.d))
  
  random.covariate <- list(l.YL = c(data.long$id.2+200, rep(NA,n2+n3+n4)), 
                           l.ZL = c(data.long$famid.2+20, rep(NA,n2+n3+n4)), 
                           z.YL = c(rep(NA,n1), data.long.z$id.2+200, rep(NA,n3+n4)), 
                           z.ZL = c(rep(NA,n1), data.long.z$famid.2+20, rep(NA,n3+n4)), 
                           r.YR = c(rep(NA,n1+n2), data.rec$id.2, rep(NA,n4)), 
                           r.ZR = c(rep(NA,n1+n2), data.rec$famid.2, rep(NA,n4)), 
                           d.YR = c(rep(NA,n1+n2+n3), data.term$id.2), 
                           d.ZR = c(rep(NA,n1+n2+n3), data.term$famid.2), 
                           d.YL = c(rep(NA,n1+n2+n3), data.term$id.2+200), 
                           d.ZL = c(rep(NA,n1+n2+n3), data.term$famid.2+20))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  # Trivariate joint model
  formula1 = Y ~ -1 + offset(c(rep(chi0,n1),rep(theta0,n2),rep(-log(b.R),n3),rep(-log(b.D),n4))) +
    l.x1 + l.x2 + offset(l.Time) +
    z.x1 + z.x2 + offset(z.Time) +
    r.x1 + r.x2 +
    d.x1 + d.x2 +
    f(d.YR, model="iidkd", order=2, n=2*200, hyper = list(theta1 = list(initial=log(L.Y[1,1]),fixed=TRUE),
                                                         theta2 = list(initial=log(L.Y[2,2]),fixed=TRUE),
                                                         theta3 = list(initial=L.Y[2,1],fixed=TRUE))) +
    f(d.YL, copy="d.YR") +
    f(r.YR, copy="d.YR", hyper = list(beta = list(initial=alpha, fixed=TRUE))) +
    f(l.YL, copy="d.YR", hyper = list(beta = list(initial=omega.l, fixed=TRUE)), n = length(sort(unique(l.YL)))) +
    f(z.YL, copy="d.YR", hyper = list(beta = list(initial=omega.z, fixed=TRUE)), n = length(sort(unique(z.YL)))) +
    f(r.ZR, model="iidkd", order=2, n=2*20, hyper = list(theta1 = list(initial=log(L.Z[1,1]),fixed=TRUE),
                                                        theta2 = list(initial=log(L.Z[2,2]),fixed=TRUE),
                                                        theta3 = list(initial=L.Z[2,1],fixed=TRUE))) +
    f(d.ZL, copy="r.ZR") +
    f(l.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.l, fixed=TRUE)), n = length(sort(unique(l.ZL)))) +
    f(z.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.z, fixed=TRUE)), n = length(sort(unique(z.ZL)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(initial=xi, fixed=TRUE))) 
  
  model3 <- INLA::inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                       data = data.f, control.compute=list(dic=FALSE,cpo=FALSE,waic=FALSE, return.marginals.predictor = FALSE,config = FALSE),
                       control.family = list(list(),
                                             list(hyper=list(size=list(initial=log(dispersion),fixed=TRUE),
                                                             prob=list(initial=-10,fixed = TRUE))),
                                             list(variant=1,hyper=list(alpha=list(initial=log(a.R),fixed=TRUE))),
                                             list(variant=1,hyper=list(alpha=list(initial=log(a.D),fixed=TRUE)))),
                       control.inla = list(int.strategy = "eb"), inla.mode = "compact", safe = TRUE, 
                       control.fixed = list(mean=list(l.x1=chi1,l.x2=chi2,
                                                      z.x1=theta1,z.x2=theta2,
                                                      r.x1=gamma1,r.x2=gamma2,
                                                      d.x1=beta1,d.x2=beta2),
                                            prec=list(l.x1=1e10,l.x2=1e10,
                                                      z.x1=1e10,z.x2=1e10,
                                                      r.x1=1e10,r.x2=1e10,
                                                      d.x1=1e10,d.x2=1e10)),
                       control.predictor = list(link=1))
  
  # model 1
  eta.c3<-model3$summary.linear.predictor$mean[n1+n2+n3+c(1:n4)]
  
  return(c(eta.c3))
}

# for t = 0
t2<-0 

data.term<-data[data$d==1,]
ind<-data.term$id.2[data.term$surv.d>t2]
n2<-length(ind)

if (n2 > 0){
  result.ind<-result.ind2<-matrix(numeric(0),nrow=200,ncol=5)
  result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  # parameter
  param<-as.numeric(est.total)
  a.D<-param[15]
  b.D<-param[12]
  
  # model 1 & 2
  for (f in 1:20){
    result<-eta.te12(est.total,1+(f-1)*10,t2,data)
    s <- 1:5
    
    # model 1
    eta.c1<-result[1,]
    result.ind[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c1,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
    
    # model 2
    eta.c2<-result[2,]
    result.ind2[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c2,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
  }
  
  write.csv(cbind(ind,result.ind[ind,]),"te-m1-t0.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2[ind,]),"te-m2-t0.csv",row.names=FALSE)
  
  # model 3
  for (i2 in 1:n2){
    i<-ind[i2]
    
    # model 3
    result3<-eta.te3(est.total,i,t2,data)
    eta.c3<-result3[1]
    result.ind3[i2,s]<-1-exp(-((t2+s)/b.D)^a.D*exp(eta.c3) + (t2/b.D)^a.D*exp(eta.c3))
    
    write.csv(cbind(ind,result.ind3),"te-m3-t0.csv",row.names=FALSE)
  }
} else {
  result.ind<-result.ind2<-result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  write.csv(cbind(ind,result.ind),"te-m1-t0.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2),"te-m2-t0.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind3),"te-m3-t0.csv",row.names=FALSE)
}

# for t = 2
t2<-2

data.term<-data[data$d==1,]
ind<-data.term$id.2[data.term$surv.d>t2]
n2<-length(ind)

if (n2 > 0){
  result.ind<-result.ind2<-matrix(numeric(0),nrow=200,ncol=5)
  result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  # parameter
  param<-as.numeric(est.total)
  a.D<-param[15]
  b.D<-param[12]
  
  # model 1 & 2
  for (f in 1:20){
    result<-eta.te12(est.total,1+(f-1)*10,t2,data)
    s <- 1:5
    
    # model 1
    eta.c1<-result[1,]
    result.ind[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c1,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
    
    # model 2
    eta.c2<-result[2,]
    result.ind2[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c2,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
  }
  
  write.csv(cbind(ind,result.ind[ind,]),"te-m1-t2.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2[ind,]),"te-m2-t2.csv",row.names=FALSE)
  
  # model 3
  for (i2 in 1:n2){
    i<-ind[i2]
    
    # model 3
    result3<-eta.te3(est.total,i,t2,data)
    eta.c3<-result3[1]
    result.ind3[i2,s]<-1-exp(-((t2+s)/b.D)^a.D*exp(eta.c3) + (t2/b.D)^a.D*exp(eta.c3))
    
    write.csv(cbind(ind,result.ind3),"te-m3-t2.csv",row.names=FALSE)
  }
} else {
  result.ind<-result.ind2<-result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  write.csv(cbind(ind,result.ind),"te-m1-t2.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2),"te-m2-t2.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind3),"te-m3-t2.csv",row.names=FALSE)
}

# for t = 5
t2<-5

data.term<-data[data$d==1,]
ind<-data.term$id.2[data.term$surv.d>t2]
n2<-length(ind)

if (n2 > 0){
  result.ind<-result.ind2<-matrix(numeric(0),nrow=200,ncol=5)
  result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  # parameter
  param<-as.numeric(est.total)
  a.D<-param[15]
  b.D<-param[12]
  
  # model 1 & 2
  for (f in 1:20){
    result<-eta.te12(est.total,1+(f-1)*10,t2,data)
    s <- 1:5
    
    # model 1
    eta.c1<-result[1,]
    result.ind[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c1,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
    
    # model 2
    eta.c2<-result[2,]
    result.ind2[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c2,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
  }
  
  write.csv(cbind(ind,result.ind[ind,]),"te-m1-t5.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2[ind,]),"te-m2-t5.csv",row.names=FALSE)
  
  # model 3
  for (i2 in 1:n2){
    i<-ind[i2]
    
    # model 3
    result3<-eta.te3(est.total,i,t2,data)
    eta.c3<-result3[1]
    result.ind3[i2,s]<-1-exp(-((t2+s)/b.D)^a.D*exp(eta.c3) + (t2/b.D)^a.D*exp(eta.c3))
    
    write.csv(cbind(ind,result.ind3),"te-m3-t5.csv",row.names=FALSE)
  }
} else {
  result.ind<-result.ind2<-result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  write.csv(cbind(ind,result.ind),"te-m1-t5.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2),"te-m2-t5.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind3),"te-m3-t5.csv",row.names=FALSE)
}

# for t = 10
t2<-10

data.term<-data[data$d==1,]
ind<-data.term$id.2[data.term$surv.d>t2]
n2<-length(ind)

if (n2 > 0){
  result.ind<-result.ind2<-matrix(numeric(0),nrow=200,ncol=5)
  result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  # parameter
  param<-as.numeric(est.total)
  a.D<-param[15]
  b.D<-param[12]
  
  # model 1 & 2
  for (f in 1:20){
    result<-eta.te12(est.total,1+(f-1)*10,t2,data)
    s <- 1:5
    
    # model 1
    eta.c1<-result[1,]
    result.ind[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c1,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
    
    # model 2
    eta.c2<-result[2,]
    result.ind2[c((1+(f-1)*10):(10+(f-1)*10)),]<-t(sapply(eta.c2,function(x) {1-exp(-((t2+s)/b.D)^a.D*exp(x) + (t2/b.D)^a.D*exp(x))}))
  }
  
  write.csv(cbind(ind,result.ind[ind,]),"te-m1-t10.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2[ind,]),"te-m2-t10.csv",row.names=FALSE)
  
  # model 3
  for (i2 in 1:n2){
    i<-ind[i2]
    
    # model 3
    result3<-eta.te3(est.total,i,t2,data)
    eta.c3<-result3[1]
    result.ind3[i2,s]<-1-exp(-((t2+s)/b.D)^a.D*exp(eta.c3) + (t2/b.D)^a.D*exp(eta.c3))
    
    write.csv(cbind(ind,result.ind3),"te-m3-t10.csv",row.names=FALSE)
  }
} else {
  result.ind<-result.ind2<-result.ind3<-matrix(numeric(0),nrow=n2,ncol=5)
  
  write.csv(cbind(ind,result.ind),"te-m1-t10.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind2),"te-m2-t10.csv",row.names=FALSE)
  write.csv(cbind(ind,result.ind3),"te-m3-t10.csv",row.names=FALSE)
}

######### Generating dynamic predictions for longitudinal outcomes #############
eta.lo12<-function(para1,i1,t1,s1){
  param<-as.numeric(para1)
  
  # longitudinal outcomes
  # logistic component
  chi0<-param[1]
  chi1<-param[2]
  chi2<-param[3]
  omega.l<-param[23]
  tau.l<-param[25]
  dispersion<-param[13]
  
  # negative binomial component
  theta0<-param[4]
  theta1<-param[5]
  theta2<-param[6]
  omega.z<-param[24]
  tau.z<-param[26]
  
  # recurrent events
  gamma1<-param[7]
  gamma2<-param[8]
  a.R<-param[14]
  b.R<-param[11]
  alpha<-param[22]
  
  # terminal event
  beta1<-param[9]
  beta2<-param[10]
  a.D<-param[15]
  b.D<-param[12]
  xi<-param[27]
  
  # random effects
  var.Y.R<-param[16]^2
  var.Y.L<-param[18]^2
  rho.Y<-param[17]
  cov.Y<-rho.Y*sqrt(var.Y.R)*sqrt(var.Y.L)
  vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),ncol=2)
  W.Y<-solve(vcov.Y)
  R.Y<-chol(W.Y)
  L.Y<-t(R.Y)
  var.Z.R<-param[19]^2
  var.Z.L<-param[21]^2
  rho.Z<-param[20]
  cov.Z<-rho.Z*sqrt(var.Z.R)*sqrt(var.Z.L)
  vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),ncol=2)
  W.Z<-solve(vcov.Z)
  R.Z<-chol(W.Z)
  L.Z<-t(R.Z)
  
  # family members' information
  famid_select <- data$famid.2[data$id.2==i1][1]
  fam_hist <- data[data$famid.2==famid_select, ]
  fam_hist$delta.R <- ifelse(fam_hist$cumgap <= t1,1,0)
  jstar.i <- aggregate(fam_hist$delta.R, by=list(fam_hist$id.2), sum)[,2]+1
  n.i <- fam_hist$n[fam_hist$time==2]
  fam_hist$jstar <- rep(jstar.i, n.i)
  
  # preparation for terminal event 
  data.term<-fam_hist[fam_hist$d==1, ]
  X.D.f<-with(data.term,cbind(x1.d,x2.d))
  id.D.f<-data.term$id.2
  famid.D.f<-data.term$famid.2
  Delta.D.f<-data.term$status.d
  T.D.f<-data.term$surv.d
  T.fh.D <- pmin(T.D.f,t1)
  delta.fh.D <- ifelse(T.D.f <= t1, Delta.D.f, 0)
  
  # preparation for recurrent events
  data.rec1<-fam_hist[fam_hist$jstar==1 & fam_hist$time==2,]
  X.R.f1<-with(data.rec1,cbind(x1.r,x2.r))
  id.R.f1<-data.rec1$id.2
  famid.R.f1<-data.rec1$famid.2
  Delta.R.f1<-data.rec1$delta.R
  T.R.f1<-rep(t1,nrow(data.rec1))
  data.rec2<-fam_hist[fam_hist$jstar>1 & fam_hist$time <= fam_hist$jstar+1,] 
  X.R.f2 <- with(data.rec2,cbind(x1.r,x2.r))
  id.R.f2<-data.rec2$id.2
  famid.R.f2<-data.rec2$famid.2
  Delta.R.f2<-data.rec2$delta.R*data.rec2$visit
  T.R.f2<-with(data.rec2, ifelse(delta.R==1, T.R, t1 - cumgap + T.R))
  T.R.f2[T.R.f2<1e-10]<-0
  data.rec<-rbind(data.rec1,data.rec2)
  X.R.f<-rbind(X.R.f1,X.R.f2)
  id.R.f<-c(id.R.f1,id.R.f2)
  famid.R.f<-c(famid.R.f1,famid.R.f2)
  Delta.R.f<-c(Delta.R.f1,Delta.R.f2)
  T.R.f<-c(T.R.f1,T.R.f2)
  
  # preparation for longitudinal outcomes
  # logistic component
  data.long<-data.rec2[data.rec2$delta.R==1,]
  X.L.f<-with(data.long,cbind(x1.l,x2.l))
  z<-rep(0,nrow(data.long))
  z[data.long$y==0]<-1
  
  # negative binomial component
  data.long.z<-data.long[data.long$y != 0,]
  X.L.f.z<-with(data.long.z,cbind(x1.l,x2.l))
  
  # other preparations
  n1<-nrow(data.long)
  n2<-nrow(data.long.z)
  n3<-nrow(data.rec)
  n4<-nrow(data.term)
  
  # construct model
  l.long<-c(z,rep(NA,n2+n3+n4+2))
  z.long <- c(rep(NA,n1),data.long.z$y,rep(NA,n3+n4+2))
  y.recu<-inla.surv(time=c(rep(NA,n1+n2),T.R.f,rep(NA,n4+2)), event=c(rep(NA,n1+n2),Delta.R.f,rep(NA,n4+2)))
  y.term<-inla.surv(time=c(rep(NA,n1+n2+n3),T.fh.D,rep(NA,2)), event=c(rep(NA,n1+n2+n3),delta.fh.D,rep(NA,2)))
  y.joint <- list(l.long, z.long, y.recu, y.term)
  
  linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n2),rep(1,n3),rep(2,n4),rep(NA,2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.Time = c(log(data.long$T.R), rep(NA,n2+n3+n4+2)), # related to logit(pi): logistic component for longitudinal outcomes 
                                 l.x1 = c(data.long$x1.l, rep(NA,n2+n3+n4),data$x1.l[data$id.2==i1][1],rep(NA,1)),
                                 l.x2 = c(data.long$x2.l, rep(NA,n2+n3+n4),data$x2.l[data$id.2==i1][1],rep(NA,1)),
                                 z.Time = c(rep(NA,n1), log(data.long.z$T.R), rep(NA,n3+n4+2)), # related to log(mu): negative binomial component for longitudinal outcomes
                                 z.x1 = c(rep(NA,n1), data.long.z$x1.l, rep(NA,n3+n4+1),data$x1.l[data$id.2==i1][1]),
                                 z.x2 = c(rep(NA,n1), data.long.z$x2.l, rep(NA,n3+n4+1),data$x2.l[data$id.2==i1][1]),
                                 r.x1 = c(rep(NA,n1+n2), data.rec$x1.r, rep(NA,n4+2)), # relative to hazard function for recurrent events
                                 r.x2 = c(rep(NA,n1+n2), data.rec$x2.r, rep(NA,n4+2)),
                                 d.x1 = c(rep(NA,n1+n2+n3), data.term$x1.d,rep(NA,2)), # related to hazard function for terminal event
                                 d.x2 = c(rep(NA,n1+n2+n3), data.term$x2.d,rep(NA,2)))
  
  random.covariate <- list(l.YL = c(data.long$id.2+200, rep(NA,n2+n3+n4),i1+200,rep(NA,1)), 
                           l.ZL = c(data.long$famid.2+20, rep(NA,n2+n3+n4),data$famid.2[data$id.2==i1][1]+20,rep(NA,1)), 
                           z.YL = c(rep(NA,n1), data.long.z$id.2+200, rep(NA,n3+n4+1),i1+200), 
                           z.ZL = c(rep(NA,n1), data.long.z$famid.2+20, rep(NA,n3+n4+1),data$famid.2[data$id.2==i1][1]+20), 
                           r.YR = c(rep(NA,n1+n2), data.rec$id.2, rep(NA,n4+2)), 
                           r.ZR = c(rep(NA,n1+n2), data.rec$famid.2, rep(NA,n4+2)), 
                           d.YR = c(rep(NA,n1+n2+n3), data.term$id.2,rep(NA,2)), 
                           d.ZR = c(rep(NA,n1+n2+n3), data.term$famid.2,rep(NA,2)), 
                           d.YL = c(rep(NA,n1+n2+n3), data.term$id.2+200,rep(NA,2)), 
                           d.ZL = c(rep(NA,n1+n2+n3), data.term$famid.2+20,rep(NA,2)))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  # Trivariate joint model
  formula1 = Y ~ -1 + offset(c(rep(chi0,n1),rep(theta0,n2),rep(-log(b.R),n3),rep(-log(b.D),n4),chi0,theta0)) +
    l.x1 + l.x2 + offset(l.Time) +
    z.x1 + z.x2 + offset(z.Time) +
    r.x1 + r.x2 +
    d.x1 + d.x2 +
    f(d.YR, model="iidkd", order=2, n=2*200, hyper = list(theta1 = list(initial=log(L.Y[1,1]),fixed=TRUE),
                                                         theta2 = list(initial=log(L.Y[2,2]),fixed=TRUE),
                                                         theta3 = list(initial=L.Y[2,1],fixed=TRUE))) +
    f(d.YL, copy="d.YR") +
    f(r.YR, copy="d.YR", hyper = list(beta = list(initial=alpha, fixed=TRUE))) +
    f(l.YL, copy="d.YR", hyper = list(beta = list(initial=omega.l, fixed=TRUE)), n = length(sort(unique(l.YL)))) +
    f(z.YL, copy="d.YR", hyper = list(beta = list(initial=omega.z, fixed=TRUE)), n = length(sort(unique(z.YL)))) +
    f(r.ZR, model="iidkd", order=2, n=2*20, hyper = list(theta1 = list(initial=log(L.Z[1,1]),fixed=TRUE),
                                                        theta2 = list(initial=log(L.Z[2,2]),fixed=TRUE),
                                                        theta3 = list(initial=L.Z[2,1],fixed=TRUE))) +
    f(d.ZL, copy="r.ZR") +
    f(l.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.l, fixed=TRUE)), n = length(sort(unique(l.ZL)))) +
    f(z.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.z, fixed=TRUE)), n = length(sort(unique(z.ZL)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(initial=xi, fixed=TRUE))) 
  
  model1 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                 data = data.f, control.compute=list(dic=FALSE,cpo=FALSE,waic=FALSE, return.marginals.predictor = FALSE,config = FALSE),
                 control.family = list(list(),
                                       list(hyper=list(size=list(initial=log(dispersion),fixed=TRUE),
                                                       prob=list(initial=-10,fixed = TRUE))),
                                       list(variant=1,hyper=list(alpha=list(initial=log(a.R),fixed=TRUE))),
                                       list(variant=1,hyper=list(alpha=list(initial=log(a.D),fixed=TRUE)))),
                 control.inla = list(int.strategy = "eb"), inla.mode = "compact", safe = TRUE, 
                 control.fixed = list(mean=list(l.x1=chi1,l.x2=chi2,
                                                z.x1=theta1,z.x2=theta2,
                                                r.x1=gamma1,r.x2=gamma2,
                                                d.x1=beta1,d.x2=beta2),
                                      prec=list(l.x1=1e10,l.x2=1e10,
                                                z.x1=1e10,z.x2=1e10,
                                                r.x1=1e10,r.x2=1e10,
                                                d.x1=1e10,d.x2=1e10)),
                 control.predictor = list(link=c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4),1,2)))
  
  # model 1
  eta.lc1<-model1$summary.linear.predictor$mean[n1+n2+n3+n4+1]+log(s1)
  pi.c1<-exp(eta.lc1)/(1+exp(eta.lc1))
  eta.zc1<-model1$summary.linear.predictor$mean[n1+n2+n3+n4+2]+log(s1)
  mu.c1<-exp(eta.zc1)
  dy.c1<-(1-pi.c1)*mu.c1
  
  # random effects
  l.ZL.1<-model1$summary.random$r.ZR$mean[famid_select+20]
  
  # model 2
  eta.lc2<-model1$summary.linear.predictor$mean[n1+n2+n3+n4+1]+log(s1)-tau.l*l.ZL.1
  pi.c2<-exp(eta.lc2)/(1+exp(eta.lc2))
  eta.zc2<-model1$summary.linear.predictor$mean[n1+n2+n3+n4+2]+log(s1)-tau.z*l.ZL.1
  mu.c2<-exp(eta.zc2)
  dy.c2<-(1-pi.c2)*mu.c2
  
  return(rbind(dy.c1,dy.c2))
}

eta.lo3<-function(para1,i1,t1,s1){
  param<-as.numeric(para1)
  
  # longitudinal outcomes
  # logistic component
  chi0<-param[1]
  chi1<-param[2]
  chi2<-param[3]
  omega.l<-param[23]
  tau.l<-param[25]
  dispersion<-param[13]
  
  # negative binomial component
  theta0<-param[4]
  theta1<-param[5]
  theta2<-param[6]
  omega.z<-param[24]
  tau.z<-param[26]
  
  # recurrent events
  gamma1<-param[7]
  gamma2<-param[8]
  a.R<-param[14]
  b.R<-param[11]
  alpha<-param[22]
  
  # terminal event
  beta1<-param[9]
  beta2<-param[10]
  a.D<-param[15]
  b.D<-param[12]
  xi<-param[27]
  
  # random effects
  var.Y.R<-param[16]^2
  var.Y.L<-param[18]^2
  rho.Y<-param[17]
  cov.Y<-rho.Y*sqrt(var.Y.R)*sqrt(var.Y.L)
  vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),ncol=2)
  W.Y<-solve(vcov.Y)
  R.Y<-chol(W.Y)
  L.Y<-t(R.Y)
  var.Z.R<-param[19]^2
  var.Z.L<-param[21]^2
  rho.Z<-param[20]
  cov.Z<-rho.Z*sqrt(var.Z.R)*sqrt(var.Z.L)
  vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),ncol=2)
  W.Z<-solve(vcov.Z)
  R.Z<-chol(W.Z)
  L.Z<-t(R.Z)
  
  # family members' information
  famid_select <- data$famid.2[data$id.2==i1][1]
  fam_hist <- data[data$famid.2==famid_select & data$id.2==i1, ]
  fam_hist$delta.R <- ifelse(fam_hist$cumgap <= t1,1,0)
  jstar.i <- aggregate(fam_hist$delta.R, by=list(fam_hist$id.2), sum)[,2]+1
  n.i <- fam_hist$n[fam_hist$time==2]
  fam_hist$jstar <- rep(jstar.i, n.i)
  
  # preparation for terminal event 
  data.term<-fam_hist[fam_hist$d==1, ]
  X.D.f<-with(data.term,cbind(x1.d,x2.d))
  id.D.f<-data.term$id.2
  famid.D.f<-data.term$famid.2
  Delta.D.f<-data.term$status.d
  T.D.f<-data.term$surv.d
  T.fh.D <- pmin(T.D.f,t1)
  delta.fh.D <- ifelse(T.D.f <= t1, Delta.D.f, 0)
  
  # preparation for recurrent events
  data.rec1<-fam_hist[fam_hist$jstar==1 & fam_hist$time==2,]
  X.R.f1<-with(data.rec1,cbind(x1.r,x2.r))
  id.R.f1<-data.rec1$id.2
  famid.R.f1<-data.rec1$famid.2
  Delta.R.f1<-data.rec1$delta.R
  T.R.f1<-rep(t1,nrow(data.rec1))
  data.rec2<-fam_hist[fam_hist$jstar>1 & fam_hist$time <= fam_hist$jstar+1,] 
  X.R.f2 <- with(data.rec2,cbind(x1.r,x2.r))
  id.R.f2<-data.rec2$id.2
  famid.R.f2<-data.rec2$famid.2
  Delta.R.f2<-data.rec2$delta.R*data.rec2$visit
  T.R.f2<-with(data.rec2, ifelse(delta.R==1, T.R, t1 - cumgap + T.R))
  T.R.f2[T.R.f2<1e-10]<-0
  data.rec<-rbind(data.rec1,data.rec2)
  X.R.f<-rbind(X.R.f1,X.R.f2)
  id.R.f<-c(id.R.f1,id.R.f2)
  famid.R.f<-c(famid.R.f1,famid.R.f2)
  Delta.R.f<-c(Delta.R.f1,Delta.R.f2)
  T.R.f<-c(T.R.f1,T.R.f2)
  
  # preparation for longitudinal outcomes
  # logistic component
  data.long<-data.rec2[data.rec2$delta.R==1,]
  X.L.f<-with(data.long,cbind(x1.l,x2.l))
  z<-rep(0,nrow(data.long))
  z[data.long$y==0]<-1
  
  # negative binomial component
  data.long.z<-data.long[data.long$y != 0,]
  X.L.f.z<-with(data.long.z,cbind(x1.l,x2.l))
  
  # other preparations
  n1<-nrow(data.long)
  n2<-nrow(data.long.z)
  n3<-nrow(data.rec)
  n4<-nrow(data.term)
  
  # construct model
  l.long<-c(z,rep(NA,n2+n3+n4+2))
  z.long <- c(rep(NA,n1),data.long.z$y,rep(NA,n3+n4+2))
  y.recu<-inla.surv(time=c(rep(NA,n1+n2),T.R.f,rep(NA,n4+2)), event=c(rep(NA,n1+n2),Delta.R.f,rep(NA,n4+2)))
  y.term<-inla.surv(time=c(rep(NA,n1+n2+n3),T.fh.D,rep(NA,2)), event=c(rep(NA,n1+n2+n3),delta.fh.D,rep(NA,2)))
  y.joint <- list(l.long, z.long, y.recu, y.term)
  
  linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n2),rep(1,n3),rep(2,n4),rep(NA,2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.Time = c(log(data.long$T.R), rep(NA,n2+n3+n4+2)), # related to logit(pi): logistic component for longitudinal outcomes 
                                 l.x1 = c(data.long$x1.l, rep(NA,n2+n3+n4),data$x1.l[data$id.2==i1][1],rep(NA,1)),
                                 l.x2 = c(data.long$x2.l, rep(NA,n2+n3+n4),data$x2.l[data$id.2==i1][1],rep(NA,1)),
                                 z.Time = c(rep(NA,n1), log(data.long.z$T.R), rep(NA,n3+n4+2)), # related to log(mu): negative binomial component for longitudinal outcomes
                                 z.x1 = c(rep(NA,n1), data.long.z$x1.l, rep(NA,n3+n4+1),data$x1.l[data$id.2==i1][1]),
                                 z.x2 = c(rep(NA,n1), data.long.z$x2.l, rep(NA,n3+n4+1),data$x2.l[data$id.2==i1][1]),
                                 r.x1 = c(rep(NA,n1+n2), data.rec$x1.r, rep(NA,n4+2)), # relative to hazard function for recurrent events
                                 r.x2 = c(rep(NA,n1+n2), data.rec$x2.r, rep(NA,n4+2)),
                                 d.x1 = c(rep(NA,n1+n2+n3), data.term$x1.d,rep(NA,2)), # related to hazard function for terminal event
                                 d.x2 = c(rep(NA,n1+n2+n3), data.term$x2.d,rep(NA,2)))
  
  random.covariate <- list(l.YL = c(data.long$id.2+200, rep(NA,n2+n3+n4),i1+200,rep(NA,1)), 
                           l.ZL = c(data.long$famid.2+20, rep(NA,n2+n3+n4),data$famid.2[data$id.2==i1][1]+20,rep(NA,1)), 
                           z.YL = c(rep(NA,n1), data.long.z$id.2+200, rep(NA,n3+n4+1),i1+200), 
                           z.ZL = c(rep(NA,n1), data.long.z$famid.2+20, rep(NA,n3+n4+1),data$famid.2[data$id.2==i1][1]+20), 
                           r.YR = c(rep(NA,n1+n2), data.rec$id.2, rep(NA,n4+2)), 
                           r.ZR = c(rep(NA,n1+n2), data.rec$famid.2, rep(NA,n4+2)), 
                           d.YR = c(rep(NA,n1+n2+n3), data.term$id.2,rep(NA,2)), 
                           d.ZR = c(rep(NA,n1+n2+n3), data.term$famid.2,rep(NA,2)), 
                           d.YL = c(rep(NA,n1+n2+n3), data.term$id.2+200,rep(NA,2)), 
                           d.ZL = c(rep(NA,n1+n2+n3), data.term$famid.2+20,rep(NA,2)))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  # Trivariate joint model
  formula1 = Y ~ -1 + offset(c(rep(chi0,n1),rep(theta0,n2),rep(-log(b.R),n3),rep(-log(b.D),n4),chi0,theta0)) +
    l.x1 + l.x2 + offset(l.Time) +
    z.x1 + z.x2 + offset(z.Time) +
    r.x1 + r.x2 +
    d.x1 + d.x2 +
    f(d.YR, model="iidkd", order=2, n=2*200, hyper = list(theta1 = list(initial=log(L.Y[1,1]),fixed=TRUE),
                                                         theta2 = list(initial=log(L.Y[2,2]),fixed=TRUE),
                                                         theta3 = list(initial=L.Y[2,1],fixed=TRUE))) +
    f(d.YL, copy="d.YR") +
    f(r.YR, copy="d.YR", hyper = list(beta = list(initial=alpha, fixed=TRUE))) +
    f(l.YL, copy="d.YR", hyper = list(beta = list(initial=omega.l, fixed=TRUE)), n = length(sort(unique(l.YL)))) +
    f(z.YL, copy="d.YR", hyper = list(beta = list(initial=omega.z, fixed=TRUE)), n = length(sort(unique(z.YL)))) +
    f(r.ZR, model="iidkd", order=2, n=2*20, hyper = list(theta1 = list(initial=log(L.Z[1,1]),fixed=TRUE),
                                                        theta2 = list(initial=log(L.Z[2,2]),fixed=TRUE),
                                                        theta3 = list(initial=L.Z[2,1],fixed=TRUE))) +
    f(d.ZL, copy="r.ZR") +
    f(l.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.l, fixed=TRUE)), n = length(sort(unique(l.ZL)))) +
    f(z.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.z, fixed=TRUE)), n = length(sort(unique(z.ZL)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(initial=xi, fixed=TRUE))) 
  
  model3 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                 data = data.f, control.compute=list(dic=FALSE,cpo=FALSE,waic=FALSE, return.marginals.predictor = FALSE,config = FALSE),
                 control.family = list(list(),
                                       list(hyper=list(size=list(initial=log(dispersion),fixed=TRUE),
                                                       prob=list(initial=-10,fixed = TRUE))),
                                       list(variant=1,hyper=list(alpha=list(initial=log(a.R),fixed=TRUE))),
                                       list(variant=1,hyper=list(alpha=list(initial=log(a.D),fixed=TRUE)))),
                 control.inla = list(int.strategy = "eb"), inla.mode = "compact", safe = TRUE, 
                 control.fixed = list(mean=list(l.x1=chi1,l.x2=chi2,
                                                z.x1=theta1,z.x2=theta2,
                                                r.x1=gamma1,r.x2=gamma2,
                                                d.x1=beta1,d.x2=beta2),
                                      prec=list(l.x1=1e10,l.x2=1e10,
                                                z.x1=1e10,z.x2=1e10,
                                                r.x1=1e10,r.x2=1e10,
                                                d.x1=1e10,d.x2=1e10)),
                 control.predictor = list(link=c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4),1,2)))
  
  # model 3
  eta.lc3<-model3$summary.linear.predictor$mean[n1+n2+n3+n4+1]+log(s1)
  pi.c3<-exp(eta.lc3)/(1+exp(eta.lc3))
  eta.zc3<-model3$summary.linear.predictor$mean[n1+n2+n3+n4+2]+log(s1)
  mu.c3<-exp(eta.zc3)
  dy.c3<-(1-pi.c3)*mu.c3
  
  return(c(dy.c3))
}

# for all t
data.long2<-data[data$d==0,]
n1<-nrow(data.long2)
dp.i.f.1c<-dp.i.f.2c<-dp.i.f.3c<-cbind(data.long2$id.2,rep(-1,n1),data.long2$y)

if (n1 > 0){
  for (i1 in 1:n1){
    i2<-data.long2$id.2[i1]
    t2<-if (data.long2$time[i1]==2) {0} else {data.long2$cumgap[i1]-data.long2$T.R[i1]}
    s2<-data.long2$T.R[i1]
    
    # result
    result<-eta.lo12(est.total,i2,t2,s2)
    result3<-eta.lo3(est.total,i2,t2,s2)
    
    # model 1
    dp.i.f.1c[i1,2]<-result[1]
    
    # model 2
    dp.i.f.2c[i1,2]<-result[2]
    
    # model 3
    dp.i.f.3c[i1,2]<-result3
    
    write.csv(dp.i.f.1c,paste0("lo-m1.csv"),row.names=FALSE)
    write.csv(dp.i.f.2c,paste0("lo-m2.csv"),row.names=FALSE)
    write.csv(dp.i.f.3c,paste0("lo-m3.csv"),row.names=FALSE)
  }} else {
    write.csv(dp.i.f.1c,paste0("lo-m1.csv"),row.names=FALSE)
    write.csv(dp.i.f.2c,paste0("lo-m2.csv"),row.names=FALSE)
    write.csv(dp.i.f.3c,paste0("lo-m3.csv"),row.names=FALSE)
  }