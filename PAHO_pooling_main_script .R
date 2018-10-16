rm(list=ls(all=TRUE))
require(reshape)
require(reshape2)
require(abind)
library(mgcv)   
library(R2jags) 
library(splines) 

####SET INPUT PARAMETERS#############################################################################################################
countries<-c('PAHO_ar','PAHO_br', 'PAHO_co', 'PAHO_dr', 'PAHO_ec', 'PAHO_gy', 'PAHO_pr') #PAHO_mx, PAHO_hr
#countries<-c('PAHO_ar','PAHO_br', 'PAHO_co',  'PAHO_ec',  'PAHO_pr') #PAHO_mx, PAHO_hr

age_group <- '2-59m' # <2m, 2-11m, 2-23m, 2-59m, 12-23m, 24-59m
hdi_level <- 'A' # Low HDI, Med HDI, Hi  HDI, A
subnational=c(0,0,0,0,0,0,0)
max.time.points=48+1 
#####################################################################################################################################
source('PAHO_pooling_source_script.R')
output_directory<- paste0(dirname(getwd()), "/Results/",hdi_level,"_",age_group, "_subnat", max(subnational),'/')
ifelse(!dir.exists(output_directory), dir.create(output_directory), FALSE)
p=N.knots+1 #intercept and slope for each time series
q=1  #2nd level predictors of intercept and slope q=1 for intercept only

###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
model_string<-"
model{

  for(i in 1:n.countries){
    for(j in 1:N.states[i]){
#################
#Observation model
#################
      w_hat[1:ts.length[i,j], j,i] ~ dmnorm(w_true[i,j, 1:ts.length[i,j]], log_rr_prec_all[1:ts.length[i,j], 1:ts.length[i,j], j,i])
#################
#Model of 'true' time series data
#################
      w_true[i,j, 1:ts.length[i,j]] ~ dmnorm(reg_mean[i,j, 1:ts.length[i,j]], w_true_cov_inv[i,j, 1:ts.length[i,j], 1:ts.length[i,j]])
      reg_mean[i,j, 1:ts.length[i,j]]<-spl.t.std[1:ts.length[i,j],]%*%beta[i,j, 1:p]
     for(k1 in 1:ts.length[i,j]){
        for(k2 in 1:ts.length[i,j]){
           w_true_cov_inv[i,j,k1,k2]<-ifelse(k1==k2, w_true_var_inv[i,j], 0)
        }
      }
##############################################################
#Second Stage Statistical Model
##############################################################
beta[i,j, 1:p] ~ dmnorm(mu1[i,j, 1:p], Sigma_inv[i, 1:p, 1:p])
for(k in 1:p){
mu1[i,j,k] <- z[i,j, 1:q]%*%gamma[i,k, 1:q]
}
}
for(k in 1:p){
  
  ###############################################################
  #Third Stage Statistical Model
  ###############################################################
  gamma[i,k, 1:q] ~ dmnorm(mu2[i,k, 1:q], Omega_inv[k, 1:q, 1:q])
    for(l in 1:q){
     mu2[i,k,l] <- w[i, 1:m]%*%theta[k,l, 1:m]
    }
  } 
Sigma_inv[i, 1:p, 1:p] ~ dwish(I_Sigma[1:p, 1:p], (p + 1))
Sigma[i, 1:p, 1:p] <- inverse(Sigma_inv[i, 1:p, 1:p])
}

#############################################################
#Remaining Prior Distributions
#############################################################
sigma2_phi_inv <- 1/(sigma_phi*sigma_phi)
sigma_phi ~ dunif(0, 1000)

for(k in 1:p){
Omega_inv[k, 1:q, 1:q] ~ dwish(I_Omega[1:q, 1:q], (q + 1))
Omega[k, 1:q, 1:q] <- inverse(Omega_inv[k, 1:q, 1:q])
for(l in 1:q){
  for(r in 1:m){
    theta[k,l,r] ~ dnorm(0, 0.0001)
  }
  }
}

}
 "


#Model Organization
model_jags<-jags.model(textConnection(model_string),
                       data=list('n.countries' = N.countries, 
                                 'N.states' = N.states, 
                                 'w_hat' = log_rr_q_all,
                                 'log_rr_prec_all' = log_rr_prec_all,  
                                 'spl.t.std' = spl.t.std, 
                                 'ts.length' = ts.length_mat,
                                 'p'=p,
                                 'q'=q,
                                 'm'=m,
                                 'w'=w,
                                 'I_Omega'= I_Omega,
                                 'I_Sigma'=I_Sigma), n.chains=2) 

#Posterior Sampling
update(model_jags, n.iter=5000)  

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("reg_mean", "beta", "w_true"),
                                thin=10,
                                n.iter=5000)
#plot(posterior_samples, ask=TRUE)


##########################################################################################
##### test for Convergence - trace plots (specify which to test in "variable.names") #####
##########################################################################################
#par(ask=TRUE)
#dev.off()
#pdf(file = "Trace Plots.pdf")
#par(mfrow = c(3,2))
#plot(posterior_samples)
#dev.off()

##################################################################
##### model fit (DIC) - must change n.chains=2 in model_jags #####
##################################################################
# dic.samples(model_jags, 500) #with pooling josh's fix (2 runs w 500 each): deviance: -717, -715; penalty=269,266, DIC=-447, -449
 #dic.samples(model_jags, 500)

# par(mfrow=c(1,1))
# caterplot(posterior_samples)

#############################################################################
##### create posterior estimates into data frames (w.true and reg.mean) #####
#############################################################################
betas.m1<-posterior_samples[[1]][,grep("^beta.*,1]",dimnames(posterior_samples[[1]])[[2]])]
betas.m2<-posterior_samples[[1]][,grep("^beta.*,2]",dimnames(posterior_samples[[1]])[[2]])]
betas.m3<-posterior_samples[[1]][,grep("^beta.*,3]",dimnames(posterior_samples[[1]])[[2]])]
betas.m4<-posterior_samples[[1]][,grep("^beta.*,4]",dimnames(posterior_samples[[1]])[[2]])]
betas.m5<-posterior_samples[[1]][,grep("^beta.*,5]",dimnames(posterior_samples[[1]])[[2]])]
betas.m6<-posterior_samples[[1]][,grep("^beta.*,6]",dimnames(posterior_samples[[1]])[[2]])]
beta.labs<-colnames(betas.m1) # labels from matrix
beta.labs.extract<- sub(".*\\[(.*)\\].*", "\\1", beta.labs, perl=TRUE)    #extract index numbers from the labels
beta.labs.extract2<-matrix(as.numeric(unlist(strsplit(beta.labs.extract, ","))),ncol=3, byrow=TRUE) #format into a matrix

preds<-array(NA,dim=c(nrow(betas.m1), nrow(spl.t.std), ncol(betas.m1)))
dimnames(preds)[[3]]<-beta.labs.extract

for(i in 1:nrow(betas.m1)){ #by iteration
  #ds.select 
  for(j in 1:ncol(betas.m1[i,, drop=FALSE] )){
   # print(i)
    #print(j)
    spl.sub<-as.matrix(spl.t.std[1:ts.length.vec[j], ,drop=FALSE])
    preds[i,1:ts.length.vec[j],j]<-(spl.sub[,1, drop=FALSE] %*% betas.m1[i,j, drop=FALSE] 
                                    +spl.sub[,2, drop=FALSE] %*% betas.m2[i,j, drop=FALSE] 
                                    +spl.sub[,3, drop=FALSE] %*% betas.m3[i,j, drop=FALSE] 
                                    +spl.sub[,4, drop=FALSE] %*% betas.m4[i,j, drop=FALSE] 
                                    +spl.sub[,5, drop=FALSE] %*% betas.m5[i,j, drop=FALSE] )
      }
}
preds.q<-as.data.frame(t(apply(preds,c(2,3),quantile, probs=c(0.5),na.rm=TRUE)))
preds.ucl<-as.data.frame(apply(preds,c(2,3),quantile, probs=c(0.975),na.rm=TRUE))
preds.lcl<-as.data.frame(apply(preds,c(2,3),quantile, probs=c(0.025),na.rm=TRUE))
saveRDS(preds, file=paste0(output_directory, "reg_mean_with_pooling bsplines.rds"))
saveRDS(state.labels, file=paste0(output_directory,"state labels.rds"))

###Remove bias term (intercept)
preds.nobias<-array(NA,dim=c(nrow(betas.m1), nrow(spl.t.std), ncol(betas.m1)))
dimnames(preds.nobias)[[3]]<-beta.labs.extract
for(i in 1:nrow(betas.m1)){ #by iteration
  #ds.select 
  for(j in 1:ncol(betas.m1[i,, drop=FALSE] )){
    # print(i)
    #print(j)
    spl.sub<-as.matrix(spl.t.std[1:ts.length.vec[j], ,drop=FALSE])
    preds.nobias[i,1:ts.length.vec[j],j]<-(spl.sub[,2, drop=FALSE] %*% betas.m2[i,j, drop=FALSE] 
                                    +spl.sub[,3, drop=FALSE] %*% betas.m3[i,j, drop=FALSE] 
                                    +spl.sub[,4, drop=FALSE] %*% betas.m4[i,j, drop=FALSE] 
                                    +spl.sub[,5, drop=FALSE] %*% betas.m5[i,j, drop=FALSE]  )
  }
}
preds.nobias.q<-as.data.frame(t(apply(preds.nobias,c(2,3),quantile, probs=c(0.5),na.rm=TRUE)))
preds.nobias.ucl<-as.data.frame(apply(preds.nobias,c(2,3),quantile, probs=c(0.975),na.rm=TRUE))
preds.nobias.lcl<-as.data.frame(apply(preds.nobias,c(2,3),quantile, probs=c(0.025),na.rm=TRUE))

saveRDS(preds.nobias, file=paste0(output_directory,"reg_mean_with_pooling bsplines nobias.rds"))

##small multiples plot
country=beta.labs.extract2[,1]
state=beta.labs.extract2[,2]
t=beta.labs.extract2[,3]
preds.nobias.q.alt<-apply(preds.nobias,c(2,3),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
preds.nobias.q.alt<-preds.nobias.q.alt[,,order(country,state )]  #sort by country, then state
pred.labs<-dimnames(preds.nobias.q.alt)[[3]] # labels from matrix
pred.labs.extract<- sub(".*\\[(.*)\\].*", "\\1", pred.labs, perl=TRUE)    #extract index numbers from the labels
pred.labs.extract2<-matrix(as.numeric(unlist(strsplit(pred.labs.extract, ","))),ncol=3, byrow=TRUE) #format into a matrix
country2<-pred.labs.extract2[,1] 
state2<-pred.labs.extract2[,2]
country.labs2<-c('AR', 'BR','CO' ,'EC', 'DR', 'GY', 'NC', 'PR') #MX, HR
cols.plot<- c('#e41a1c','#377eb8','#4daf4a','#984ea3', '#ff7f00','#a65628','#f781bf', 'grey')
dim1<- max(1,floor(sqrt(length(countries))))
dim2<- max(1,ceiling(sqrt(length(countries))))
if(dim1+dim2<length(countries)){dim1=dim1+ 1}
tiff(paste0(output_directory,'small multiples nobias.tiff'), width = 7, height = 4, units = "in",res=200)
par(mfrow=c(dim1,dim2) , mar=c(2,2,0.6,0))
for (i in 1:dim(preds.nobias.q.alt )[3]){
  y.all<- t(exp(preds.nobias.q.alt[,,i]))
  keep.obs<-as.data.frame(complete.cases(y.all))[,1]
  y.all<-y.all[keep.obs,]
  n.times<- nrow(y.all)
  xx<- c(1:n.times, rev(1:n.times))
  yy<- c(y.all[,1], rev(y.all[,3] ))
  matplot(t(exp(preds.nobias.q.alt[,,i])), bty='l', type='l', col='white', lty=c(2,1,2), xlim=c(0,48), ylim=c(0.5,2), axes=F)
  polygon(xx, yy, col='gray90', border=NA)
  points( y.all[,2], col=cols.plot[country2[i]], type='l')
  abline(h=1, lty=2, col='black', lwd=1)
  abline(h=c(0.8,1.2), lty=2, col='gray', lwd=0.5)
  Axis(side=1, labels=FALSE, at=c(0,12,24,36,48), tck=-0.005,  col='gray')
  Axis(side=2, labels=FALSE, at=c(0.5,1,2.0), tck=-0.005, col='gray')
  if(state2[i]==1){
    title(country.labs2[country2[i]], col.main=cols.plot[country2[i]] )
    }
}
dev.off()



#########################
##PLOTS FOR FITTED VALUES, INCLUDING INTERCEPT
tiff(paste0(output_directory,'ucl and LCL smooth with pooling bsplines.tiff'), width = 7, height = 4, units = "in",res=200)
    par(mfrow=c(1,2))
    matplot(preds.ucl, type='l', bty='l',col='gray')
    title("Upper CrI")
    abline(h=0, col='red',lty=3)
    matplot(preds.lcl, type='l', bty='l',col='gray')
    title("Lower CrI")
    abline(h=0, col='red',lty=3)
dev.off()

country=beta.labs.extract2[,1]
state=beta.labs.extract2[,2]
t=beta.labs.extract2[,3]
tiff(paste0(output_directory,'all_HDI.tiff'), width = 10.5, height = 2.5, units = "in",res=200) 
    par(mfrow=c(1,N.countries))
    for(i in 1: N.countries){
      print(i)
      ds.select<-t(preds.q[country==i,,drop=FALSE])
      
      N.states<-ncol(ds.select)
      log_rr_prec_state_mean<-rep(NA, times=N.states)
      for (m in 1: N.states){
        log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
      }
      inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
      inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
      col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
      matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col=col.plot,lty=1, ylim=c(-0.5,0.5), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
      abline(h=0)
      title(countries[i])
    }
dev.off()   
      
tiff(paste0(output_directory,'ucl and LCL smooth with pooling bsplines nobias.tiff'), width = 7, height = 4, units = "in",res=200)
      par(mfrow=c(1,2))
      matplot(preds.nobias.ucl, type='l', bty='l',col='gray')
      title("Upper CrI")
      abline(h=0, col='red',lty=3)
      matplot(preds.nobias.lcl, type='l', bty='l',col='gray')
      title("Lower CrI")
      abline(h=0, col='red',lty=3)
      
dev.off()
  
tiff(paste0(output_directory,'smooth with pooling bspline nobias.tiff'), width = 10.5, height = 4, units = "in",res=200)
  country=beta.labs.extract2[,1]
  state=beta.labs.extract2[,2]
  t=beta.labs.extract2[,3]
    par(mfrow=c(1,N.countries))
  for(i in 1: N.countries){
    print(i)
    ds.select<-t(preds.nobias.q[country==i,,drop=FALSE])
    
    N.states<-ncol(ds.select)
    log_rr_prec_state_mean<-rep(NA, times=N.states)
    for (m in 1: N.states){
      log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
    }  
        inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
        inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
        col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
        matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col=col.plot,lty=1, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
        abline(h=0)
        title(countries[i])
         }
dev.off()


save(list = ls(all.names = TRUE),file=paste0(paste0(output_directory,"pooling no covars all states b splines.RData"),envir = .GlobalEnv))

