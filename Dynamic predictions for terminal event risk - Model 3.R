rm(list=ls())
library(INLA)
inla.setOption(inla.mode="experimental")
inla.setOption(num.threads=1) 

# import the data
data<-read.table("NFdat_newerdate.txt", stringsAsFactors=FALSE)

# preparation for longitudinal outcomes
# remove the observations with num >= 99
data$num[data$num>98] = NA

n_1<-nrow(data)-1
for (i in 1:n_1){
  if (!is.na(data$polyps[i]) & data$polyps[i] == 1 & is.na(data$num[i])) {
    data$num[i] = 1
  }
  if (!is.na(data$polyps[i]) & data$polyps[i] == 0 & is.na(data$num[i])) {
    data$num[i] = 0
  }
  if (is.na(data$polyps[i]) & is.na(data$num[i]) & data$id[i] == data$id[i+1] & data$polyps0[i+1] == 1) {
    data$num[i] = 1
    data$polyps[i] = 1
  }
  if (is.na(data$polyps[i]) & is.na(data$num[i]) & data$id[i] == data$id[i+1] & data$polyps0[i+1] == 0) {
    data$num[i] = 0
    data$polyps[i] = 0
  }
}

# cumulative gap at previous visit
data$cumgap0 <- 0
for(i in 1:177){
  n <- sum(data$id==i)
  data$cumgap0[data$id==i] <- c(0, data$cumgap[data$id==i][-n])
}

# import sampled parameters
param_sample<-read.csv("param_sample.csv")

# generate the linear predictor (eta) in terminal event process by using the conditional approach
eta.con<-function(para1,i1,t1){
  # import sampled parameters
  param<-as.vector(unlist(para1))
  
  # longitudinal outcomes
  # logistic component
  chi0<-param[1]
  chi1<-param[2]
  chi2<-param[3]
  chi3<-param[4]
  chi4<-param[5]
  omega.l<-param[31]
  tau.l<-param[33]
  dispersion<-param[11]
  
  # negative binomial component
  theta0<-param[6]
  theta1<-param[7]
  theta2<-param[8]
  theta3<-param[9]
  theta4<-param[10]
  omega.z<-param[32]
  tau.z<-param[34]
  
  # recurrent events
  gamma1<-param[12]
  gamma2<-param[13]
  gamma3<-param[14]
  gamma4<-param[15]
  a.R<-param[16]
  b.R<-param[17]
  alpha<-param[30]
  
  # terminal event
  beta1<-param[18]
  beta2<-param[19]
  beta3<-param[20]
  beta4<-param[21]
  a.D<-param[22]
  b.D<-param[23]
  xi<-param[35]
  
  # random effects
  var.Y.R<-param[24]^2
  var.Y.L<-param[26]^2
  rho.Y<-param[25]
  cov.Y<-rho.Y*sqrt(var.Y.R)*sqrt(var.Y.L)
  vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),ncol=2)
  W.Y<-solve(vcov.Y)
  R.Y<-chol(W.Y)
  L.Y<-t(R.Y)
  var.Z.R<-param[27]^2
  var.Z.L<-param[29]^2
  rho.Z<-param[28]
  cov.Z<-rho.Z*sqrt(var.Z.R)*sqrt(var.Z.L)
  vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),ncol=2)
  W.Z<-solve(vcov.Z)
  R.Z<-chol(W.Z)
  L.Z<-t(R.Z)
  
  # family members' information
  date_select <- data$datecx1_num[data$id==i1][1]
  famid_select <- data$famid[data$id==i1][1]
  fam_select <- data[data$famid==famid_select, ]
  fam_select$dur <- (date_select-fam_select$datecx1_num)/365 
  fam_hist<-fam_select[(fam_select$dur + t1)>=0 & fam_select$id==i1,] # without family history
  fam_hist$delta.R <- ifelse(fam_hist$cumgap <= (fam_hist$dur + t1),1,0)
  jstar.i <- aggregate(fam_hist$delta.R, by=list(fam_hist$id), sum)[,2]+1
  n.i <- fam_hist$n[fam_hist$time==2]
  fam_hist$jstar <- rep(jstar.i, n.i)
  
  # preparation for terminal event 
  data.term<-fam_hist[fam_hist$type=="crc", ]
  X.D.f<-with(data.term,cbind(sex01, adenoma00, age00_trans, pcrc1age_trans))
  id.D.f<-data.term$id
  famid.D.f<-data.term$famid
  Delta.D.f<-data.term$crc1
  T.D.f<-data.term$crcgap
  dur.f<-data.term$dur
  T.fh.D <- pmin(T.D.f,t1+dur.f)
  delta.fh.D <- ifelse(T.D.f <= t1+dur.f, Delta.D.f, 0)

  # preparation for recurrent events
  data.rec1<-fam_hist[fam_hist$jstar==1 & fam_hist$time==2,]
  X.R.f1<-with(data.rec1,cbind(sex01, adenoma0, age0_trans, pcrc1age_trans))
  id.R.f1<-data.rec1$id
  famid.R.f1<-data.rec1$famid
  Delta.R.f1<-data.rec1$delta.R
  T.R.f1<-data.rec1$dur+t1
  data.rec2<-fam_hist[fam_hist$jstar>1 & fam_hist$time <= fam_hist$jstar+1,] 
  X.R.f2 <- with(data.rec2,cbind(sex01, adenoma0, age0_trans, pcrc1age_trans))
  id.R.f2<-data.rec2$id
  famid.R.f2<-data.rec2$famid
  Delta.R.f2<-data.rec2$delta.R*data.rec2$visit
  T.R.f2<-with(data.rec2, ifelse(delta.R==1, gap, dur+t1 - cumgap + gap))
  T.R.f2[T.R.f2<1e-10]<-0
  data.rec<-rbind(data.rec1,data.rec2)
  X.R.f<-rbind(X.R.f1,X.R.f2)
  id.R.f<-c(id.R.f1,id.R.f2)
  famid.R.f<-c(famid.R.f1,famid.R.f2)
  Delta.R.f<-c(Delta.R.f1,Delta.R.f2)
  T.R.f<-c(T.R.f1,T.R.f2)

  # preparation for longitudinal outcomes
  # logistic component
  data.long<-data.rec2[data.rec2$delta.R==1 & !is.na(data.rec2$num),]
  X.L.f<-with(data.long,cbind(sex01, age00_trans, pcrc1age_trans, cumgap0))
  z<-rep(0,nrow(data.long))
  z[data.long$num==0]<-1
  
  # negative binomial component
  data.long.z<-data.long
  X.L.f.z<-with(data.long.z,cbind(sex01, age00_trans, pcrc1age_trans, cumgap0))
  
  # other preparations
  n1<-nrow(data.long)
  n2<-nrow(data.long.z)
  n3<-nrow(data.rec)
  n4<-nrow(data.term)
  n5<-length(unique(data.term$id))
  n6<-length(unique(data.term$famid))
  
  # construct model
  l.long<-c(z,rep(NA,n2+n3+n4+1))
  z.long <- c(rep(NA,n1),data.long.z$num,rep(NA,n3+n4+1))
  y.recu<-inla.surv(time=c(rep(NA,n1+n2),T.R.f/max(T.R.f),rep(NA,n4)), event=c(rep(NA,n1+n2),Delta.R.f,rep(NA,n4+1)))
  y.term<-inla.surv(time=c(rep(NA,n1+n2+n3),T.fh.D/max(T.fh.D)), event=c(rep(NA,n1+n2+n3),delta.fh.D,rep(NA,1)))
  y.joint <- list(l.long, z.long, y.recu, y.term)
  
  linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n2),rep(1,n3),rep(2,n4+1))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.Time = c(log(data.long$gap), rep(NA,n2+n3+n4+1)), # related to logit(pi): logistic component for longitudinal outcomes 
                                 l.Sex = c(data.long$sex01, rep(NA,n2+n3+n4+1)),
                                 l.Age = c(data.long$age00_trans, rep(NA,n2+n3+n4+1)),
                                 l.Pcrc1age = c(data.long$pcrc1age_trans, rep(NA,n2+n3+n4+1)),
                                 l.Cumgap = c(data.long$cumgap0, rep(NA,n2+n3+n4+1)),
                                 z.Time = c(rep(NA,n1), log(data.long.z$gap), rep(NA,n3+n4+1)), # related to log(mu): negative binomial component for longitudinal outcomes
                                 z.Sex = c(rep(NA,n1), data.long.z$sex01, rep(NA,n3+n4+1)),
                                 z.Age = c(rep(NA,n1), data.long.z$age00_trans, rep(NA,n3+n4+1)),
                                 z.Pcrc1age = c(rep(NA,n1), data.long.z$pcrc1age_trans, rep(NA,n3+n4+1)),
                                 z.Cumgap = c(rep(NA,n1), data.long.z$cumgap0, rep(NA,n3+n4+1)),   
                                 r.Sex = c(rep(NA,n1+n2), data.rec$sex01, rep(NA,n4+1)), # relative to hazard function for recurrent events
                                 r.Adenoma = c(rep(NA,n1+n2), data.rec$adenoma0, rep(NA,n4+1)),
                                 r.Age = c(rep(NA,n1+n2), data.rec$age0_trans, rep(NA,n4+1)),
                                 r.Pcrc1age = c(rep(NA,n1+n2), data.rec$pcrc1age_trans, rep(NA,n4+1)),
                                 d.Sex = c(rep(NA,n1+n2+n3), data.term$sex01,data$sex01[data$id==i1][1]), # related to hazard function for terminal event
                                 d.Adenoma = c(rep(NA,n1+n2+n3), data.term$adenoma00,data$adenoma00[data$id==i1][1]),
                                 d.Age = c(rep(NA,n1+n2+n3), data.term$age00_trans,data$age00_trans[data$id==i1][1]),
                                 d.Pcrc1age = c(rep(NA,n1+n2+n3), data.term$pcrc1age_trans,data$pcrc1age_trans[data$id==i1][1]))
  
  random.covariate <- list(l.YL = c(data.long$id+177, rep(NA,n2+n3+n4+1)), 
                           l.ZL = c(data.long$famid+18, rep(NA,n2+n3+n4+1)), 
                           z.YL = c(rep(NA,n1), data.long.z$id+177, rep(NA,n3+n4+1)), 
                           z.ZL = c(rep(NA,n1), data.long.z$famid+18, rep(NA,n3+n4+1)), 
                           r.YR = c(rep(NA,n1+n2), data.rec$id, rep(NA,n4+1)), 
                           r.ZR = c(rep(NA,n1+n2), data.rec$famid, rep(NA,n4+1)), 
                           d.YR = c(rep(NA,n1+n2+n3), data.term$id,i1), 
                           d.ZR = c(rep(NA,n1+n2+n3), data.term$famid,famid_select), 
                           d.YL = c(rep(NA,n1+n2+n3), data.term$id+177,i1+177), 
                           d.ZL = c(rep(NA,n1+n2+n3), data.term$famid+18,famid_select+18))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  # Trivariate joint model
  formula1 = Y ~ -1 + offset(c(rep(chi0,n1),rep(theta0,n2),rep(-log(b.R),n3),rep(-log(b.D),n4+1))) +
    l.Sex + l.Age + l.Pcrc1age + l.Cumgap + offset(l.Time) +
    z.Sex + z.Age + z.Pcrc1age + z.Cumgap + offset(z.Time) +
    r.Sex + r.Adenoma + r.Age + r.Pcrc1age +
    d.Sex + d.Adenoma + d.Age + d.Pcrc1age +
    f(d.YR, model="iidkd", order=2, n=2*177, hyper = list(theta1 = list(initial=log(L.Y[1,1]),fixed=TRUE),
                                                          theta2 = list(initial=log(L.Y[2,2]),fixed=TRUE),
                                                          theta3 = list(initial=L.Y[2,1],fixed=TRUE))) +
    f(d.YL, copy="d.YR") +
    f(r.YR, copy="d.YR", hyper = list(beta = list(initial=alpha, fixed=TRUE))) +
    f(l.YL, copy="d.YR", hyper = list(beta = list(initial=omega.l, fixed=TRUE)), n = length(sort(unique(l.YL)))) +
    f(z.YL, copy="d.YR", hyper = list(beta = list(initial=omega.z, fixed=TRUE)), n = length(sort(unique(z.YL)))) +
    f(r.ZR, model="iidkd", order=2, n=2*18, hyper = list(theta1 = list(initial=log(L.Z[1,1]),fixed=TRUE),
                                                         theta2 = list(initial=log(L.Z[2,2]),fixed=TRUE),
                                                         theta3 = list(initial=L.Z[2,1],fixed=TRUE))) +
    f(d.ZL, copy="r.ZR") +
    f(l.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.l, fixed=TRUE)), n = length(sort(unique(l.ZL)))) +
    f(z.ZL, copy="r.ZR", hyper = list(beta = list(initial=tau.z, fixed=TRUE)), n = length(sort(unique(z.ZL)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(initial=xi, fixed=TRUE))) 
  
  model1 <- inla(formula1, family = c("binomial","zeroinflatednbinomial1","weibullsurv","weibullsurv"),
                 data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE, return.marginals.predictor = TRUE,config = TRUE),
                 control.family = list(list(),
                                       list(hyper=list(size=list(initial=log(dispersion),fixed=TRUE),
                                                       prob=list(initial=-10,fixed = TRUE))),
                                       list(variant=1,hyper=list(alpha=list(initial=log(a.R),fixed=TRUE))),
                                       list(variant=1,hyper=list(alpha=list(initial=log(a.D),fixed=TRUE)))),
                 control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE,
                 control.fixed = list(mean=list(l.Sex=chi1,l.Age=chi2,l.Pcrc1age=chi3,l.Cumgap=chi4,
                                                z.Sex=theta1,z.Age=theta2,z.Pcrc1age=theta3,z.Cumgap=theta4,
                                                r.Sex=gamma1,r.Adenoma=gamma2,r.Age=gamma3,r.Pcrc1age=gamma4,
                                                d.Sex=beta1,d.Adenoma=beta2,d.Age=beta3,d.Pcrc1age=beta4),
                                      prec=list(l.Sex=1e10,l.Age=1e10,l.Pcrc1age=1e10,l.Cumgap=1e10,
                                                z.Sex=1e10,z.Age=1e10,z.Pcrc1age=1e10,z.Cumgap=1e10,
                                                r.Sex=1e10,r.Adenoma=1e10,r.Age=1e10,r.Pcrc1age=1e10,
                                                d.Sex=1e10,d.Adenoma=1e10,d.Age=1e10,d.Pcrc1age=1e10)),
                 control.predictor = list(link=1))
  
  # model 3
  eta.c3<-model1$summary.linear.predictor$mean[n1+n2+n3+which(id.D.f==i1)]
  
  return(c(eta.c3))
}

# for t = 0
t2<-0 
data.term<-data[data$type=="crc",]
ind<-data.term$id[data.term$crcgap>t2]
n2<-length(ind)
m2<-5
dp.i.f.3c<-matrix(-1,nrow=n2,ncol=15)
for (i2 in 1:n2){
  i<-ind[i2]
  result.ind3<-matrix(-1,nrow=m2,ncol=15)
  
  for (j in 1:m2){
    # parameter
    param<-as.vector(unlist(param_sample[j,]))
    a.D<-param[22]
    b.D<-param[23]

    # result
    result<-eta.con(param_sample[j,],i,t2)
    s <- 1:15
    
    # model 3
    eta.c3<-result[1]
    result.ind3[j,s]<-1-exp(-((t2+s)/b.D)^a.D*exp(eta.c3) + (t2/b.D)^a.D*exp(eta.c3))

    print(c(i,j))
  }
  
  dp.i.f.3c[i2,]<-apply(result.ind3,2,median,na.rm=TRUE)
  
  write.csv(cbind(ind,dp.i.f.3c),"dp-te-3c-0.csv",row.names=FALSE)
}




