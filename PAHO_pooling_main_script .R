
#Loop trough all of the age categories
#for(agegrp.select in c('<2m', '2-11m', '2-23m', '2-59m', '12-23m', '24-59m')){
for(agegrp.select in c('2-11m', '<2m', '12-23m', '24-59m','2-59m')){
  rm(list=ls()[-which(ls() %in% c('agegrp.select'))]) 
  
print(agegrp.select)
require(reshape)
require(reshape2)
require(abind)
library(mgcv)   
library(R2jags) 
library(splines) 
library(lubridate)

####SET INPUT PARAMETERS#############################################################################################################
countries<-c('PAHO_ar','PAHO_br', 'PAHO_co','PAHO_dr',  'PAHO_ec','PAHO_hr', 'PAHO_mxA','PAHO_nc','PAHO_pr') # PAHO_hr 'PAHO_gy',PAHO_nc

age_group <- agegrp.select # <2m, 2-11m, 2-23m, 2-59m, 12-23m, 24-59m
hdi_level <- 'A' # Low HDI, Med HDI, Hi  HDI, A
subnational=rep(0, length(countries))
max.time.points=48
pre.vax.time<-12 #how many to use when anchoring at t=0
tot_time<-max.time.points+pre.vax.time
#####################################################################################################################################
source('PAHO_pooling_source_script.R')
output_directory<- paste0(dirname(getwd()), "/Results CP/",hdi_level,"_",age_group, "_nat_noMVN", max(subnational),'/')
output_directory<-gsub("<2", "u2", output_directory)
ifelse(!dir.exists(output_directory), dir.create(file.path(output_directory)), FALSE)

matplot(log_rr_q_all[,1,], type='l', bty='l', col=1:N.countries, ylim=c(-0.5,0.5))
legend(x=30, y=-1, countries, col=1:N.countries, lty=2)
abline(h=0, v=pre.vax.time)

##test Mexico
 # mod1<-lm(log_rr_q_all[,1,7] ~ spl.t.std[,2]+spl.t.std[,3] +spl.t.std[,4]+spl.t.std[,5])
 # pred1<-predict(mod1)
 # plot(pred1-coef(mod1)[1], ylim=c(min(pred1),0.2))

#Check if country has data, and if not, remove it from the array
nodata.country<-is.na(log_rr_q_all[1,1,])
log_rr_q_all<-log_rr_q_all[,,!nodata.country, drop=FALSE]
log_rr_prec_all<-log_rr_prec_all[,,,!nodata.country, drop=FALSE]
log_rr_prec_point_all<-log_rr_prec_point_all[,,!nodata.country, drop=FALSE]

countries<-countries[!nodata.country]
saveRDS(countries, file=paste0(output_directory, "country_labs.rds"))

N.countries<-length(countries)
###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
time.index<-c(rep(0,times=pre.vax.time), seq.int(from=1, to=max.time.points, by=1))/max.time.points

model_string<-"
model{
for(i in 1:n.countries){
for(j in 1:n.states[i]){
for(v in 1:ts.length[i,j]){       

##################
#Observation model
##################
w_hat[v, j,i] ~ dnorm(w_true[i,j, v], log_rr_prec_all[ v, j,i])

#################################
#Model of 'true' time series data
#################################
w_true[i,j, v] ~ dnorm(reg_mean[i,j,v], w_true_prec_inv[i,j])

##################### 
#CHANGE POINT MODEL #
#####################
reg_mean[i,j,v]<-(beta[i,j,1] +
           step(time.index[v] - cp1[i,j])*(1 - step(time.index[v] - cp2[i,j]))*beta[i,j,2]*(time.index[v] - cp1[i,j]) 
          +step(time.index[v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j]))
reg_unbias[i,j,v]<-(step(time.index[v] - cp1[i,j])*(1 - step(time.index[v] - cp2[i,j]))*beta[i,j,2]*(time.index[v] - cp1[i,j]) 
              +step(time.index[v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j]))
}
#slope[i,j]<- -exp(beta[i,j,2]) #Ensures slope is negative

w_true_prec_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
w_true_sd[i,j] ~ dunif(0, 100)

cp1[i,j]<-exp(beta[i,j,3])
cp2.add[i,j]<-exp(beta[i,j,4])
cp2[i,j]<-cp1[i,j] +cp2.add[i,j]  + 1/max.time.points   #ensure Cp2 is at least 1 unit after CP1
###########################################################
#Second Stage Statistical Model
###########################################################
beta[i,j, 1:4] ~ dmnorm(lambda[1:4], Omega_inv[1:4, 1:4])
}
########################################################
#Third Stage Statistical Model
########################################################
#gamma[i, 1:4] ~ dmnorm(lambda[1:4], Sigma_inv[1:4, 1:4])
}
#######################################################
#Remaining Prior Distributions
#######################################################
Omega_inv[1:4, 1:4] ~ dwish(I_Omega[1:4, 1:4], (4 + 1))
Omega[1:4, 1:4]<-inverse(Omega_inv[1:4, 1:4])
# Sigma_inv[1:4, 1:4] ~ dwish(I_Sigma[1:4, 1:4], (4 + 1))
# Sigma[1:4, 1:4]<-inverse(Sigma_inv[1:4, 1:4])

for(j in c(1:4)){
lambda[j] ~ dnorm(0, 1e-4)
}

}
"
#Model Organization
model_jags<-jags.model(textConnection(model_string),
                       data=list('n.countries' = N.countries, 
                                 'n.states' = N.states, 
                                 'w_hat' = log_rr_q_all,
                                 'log_rr_prec_all' = log_rr_prec_point_all, 
                                 'ts.length' = ts.length_mat,
                                 'I_Omega'= I_Omega,
                                 'max.time.points'=max.time.points,
                                 'time.index'=time.index,
                                 'I_Sigma'=I_Sigma), n.chains=2, n.adapt=2000) 

#Posterior Sampling
update(model_jags, n.iter=10000)  

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("reg_mean",'reg_unbias' ,'cp1','cp2',"beta",'lambda'),
                                thin=10,
                                n.iter=10000)
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
#############################################################################
##### create posterior estimates into data frames (w.true and reg.mean) #####
#############################################################################
x.func <- function(x){ sub(".*\\[(.*)\\].*", "\\1", x, perl=TRUE)} 

beta1<-posterior_samples[[1]][,grep("^beta.*,1]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
beta1.lab<-x.func(dimnames(beta1)[[2]]) 
beta1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta1.lab.spl)<-c('country','state')

#Extract first change point time:
beta3<-posterior_samples[[1]][,grep("^beta.*,3]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
cp1<-max.time.points*exp(beta3) #Intercept
beta3.lab<-x.func(dimnames(beta3)[[2]]) 
beta3.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta3.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta3.lab.spl)<-c('country','state')
#for(i in c(1:ncol(beta3))){hist(beta3[,i])}
quant.cp1<-t(max.time.points*exp(apply(beta3,2,quantile, probs=c(0.025,0.5,0.975))))
var.cp1<-apply(beta3,2,var)
cp1.lab<-x.func(dimnames(quant.cp1)[[1]]) 
cp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(cp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(cp1.lab.spl)<-c('country','state')
par(mfrow=c(2,2))
plot(y=1:nrow(quant.cp1), x=quant.cp1[,'50%'], bty='l')
library(ggplot2)
plot.cp<-cbind.data.frame('strata'=1:nrow(quant.cp1),'median.cp'=quant.cp1[,'50%'],'lcl.cp'=quant.cp1[,'2.5%'],'ucl.cp'=quant.cp1[,'97.5%'], 'inv.var.cp'=1/var.cp1,beta3.lab.spl)
plot.cp<-plot.cp[order(plot.cp$country),]
plot.cp$order2<-1:nrow(plot.cp)
plot.cp$country2<-NA
for(i in 1:length(countries)){plot.cp$country2[plot.cp$country==i] <-countries[i] }
tiff(paste0(output_directory,'cp by country.tiff'), width = 7, height = 8, units = "in",res=200)
cp.plot<-ggplot(data=plot.cp, aes(x=median.cp, y=order2, color=country2)) +
  geom_point(aes(size=inv.var.cp)) +
  scale_size_continuous(range=c(1,15)) +
  theme_bw()+
  guides( size = FALSE)+
  # theme(legend.position = "none")+
  scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',  '#6a3d9a'))
print(cp.plot)
dev.off()
saveRDS(plot.cp, file=paste0(output_directory,'CP1 pool model.rds'))


############################
#Extract slope
beta2<-posterior_samples[[1]][,grep("^beta.*,2]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
slp1<-beta2 #Intercept
beta2.lab<-x.func(dimnames(beta2)[[2]]) 
beta2.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta2.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta2.lab.spl)<-c('country','state')
#for(i in c(1:9)){hist(beta2[,i])}
quant.slp1<-t(apply(beta2,2,quantile, probs=c(0.025,0.5,0.975)))
var.slp1<-apply(beta2,2,var)
slp1.lab<-x.func(dimnames(quant.slp1)[[1]]) 
slp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(slp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(slp1.lab.spl)<-c('country','state')
par(mfrow=c(2,2))
plot(y=1:nrow(quant.slp1), x=quant.slp1[,'50%'], bty='l')
library(ggplot2)
plot.slp<-cbind.data.frame('strata'=1:nrow(quant.slp1),'median.slp'=quant.slp1[,'50%'], 'inv.var.slp'=1/var.slp1,beta2.lab.spl)
plot.slp<-plot.slp[order(plot.slp$country),]
plot.slp$order2<-1:nrow(plot.slp)
plot.slp$country2<-NA
for(i in 1:length(countries)){plot.slp$country2[plot.slp$country==i] <-countries[i] }
tiff(paste0(output_directory,'country.slope.tiff'), width = 7, height = 8, units = "in",res=200)
slp.plot<-ggplot(data=plot.slp, aes(x=median.slp, y=order2, color=country2)) +
  geom_point(aes(size=inv.var.slp)) +
  scale_size_continuous(range=c(1,15)) +
  theme_bw()+
  guides( size = FALSE)+
  # theme(legend.position = "none")+
  scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',  '#6a3d9a'))
print(slp.plot)
dev.off()
saveRDS(plot.slp, file=paste0(output_directory,"slope pool model.rds"))

##########################################

#Relationship between change point location and slope
tiff(paste0(output_directory,'cp vs slope.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(quant.slp1[,'50%'], quant.cp1[,'50%'], col=beta2.lab.spl$country,bty='l', ylab="Change point", xlab="Slope")
dev.off()

##melt and cast predicted values into 4D array N,t,i,j array
reg_mean<-posterior_samples[[1]][,grep("reg_mean",dimnames(posterior_samples[[1]])[[2]])]
pred1<-reg_mean
pred.indices<- x.func(dimnames(pred1)[[2]]) 
pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
pred.indices.spl<-as.data.frame(pred.indices.spl)
names(pred.indices.spl)<- c('country','state','time')
reg_mean2<-cbind.data.frame(pred.indices.spl,t(reg_mean))
reg_mean_m<-melt(reg_mean2, id=c('time','country','state'))
reg_mean_c<-acast(reg_mean_m, variable~time~country~state)
preds.q<-apply(reg_mean_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.q)[[2]]<- as.numeric(as.character(dimnames(preds.q)[[2]]))

#Unbias
reg_unbias<-posterior_samples[[1]][,grep("reg_unbias",dimnames(posterior_samples[[1]])[[2]])]
pred1<-reg_unbias
pred.indices<- x.func(dimnames(pred1)[[2]]) 
pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
pred.indices.spl<-as.data.frame(pred.indices.spl)
names(pred.indices.spl)<- c('country','state','time')
reg_unbias2<-cbind.data.frame(pred.indices.spl,t(reg_unbias))
reg_unbias_m<-melt(reg_unbias2, id=c('time','country','state'))
reg_unbias_c<-acast(reg_unbias_m, variable~time~country~state)
reg_unbias_c<-reg_unbias_c[,order(as.numeric(dimnames(reg_unbias_c)[[2]])),order(as.numeric(dimnames(reg_unbias_c)[[3]])),order(as.numeric(dimnames(reg_unbias_c)[[4]]))]
preds.unbias.q<-apply(reg_unbias_c,c(2,3),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.unbias.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))

tiff(paste0(output_directory,'country.trend.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(4,2,1,1))
for(i in 1:length(countries)){
 # for(j in 1:N.states[i]){
    plot.data<-t(preds.unbias.q[,,i])
    matplot( ((1:tot_time)-pre.vax.time), plot.data,type='l',yaxt='n', xlim=c(0, max.time.points), xlab='months post-PCV introduction', ylim=c(-0.7,0.4), col='gray', lty=c(2,1,2), bty='l')
    abline(h=0)
    axis(side=2, at=c(-0.7,-0.35,0,0.35,0.7), las=1,labels=round(exp(c(-0.7,-0.35,0,0.35,0.7)),1 ))
    # abline(v=0)
    title(countries[i])
  }
dev.off()
saveRDS(preds.unbias.q, file=paste0(output_directory,"reg_mean_with_pooling CP nobias.rds"))



saveRDS(reg_unbias_c, file=paste0(output_directory, "reg_mean_unbias_with_pooling CP.rds"))
saveRDS(state.labels, file=paste0(output_directory,"state labels.rds"))
saveRDS(posterior_samples, file=paste0(output_directory, "posterior samples pooling with CP.rds"))

# 
# ##small multiples plot
# country=beta.labs.extract2[,1]
# state=beta.labs.extract2[,2]
# t=beta.labs.extract2[,3]
# preds.nobias.q.alt<-apply(preds.nobias,c(2,3),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
# preds.nobias.q.alt<-preds.nobias.q.alt[,,order(country,state )]  #sort by country, then state
# pred.labs<-dimnames(preds.nobias.q.alt)[[3]] # labels from matrix
# pred.labs.extract<- sub(".*\\[(.*)\\].*", "\\1", pred.labs, perl=TRUE)    #extract index numbers from the labels
# pred.labs.extract2<-matrix(as.numeric(unlist(strsplit(pred.labs.extract, ","))),ncol=3, byrow=TRUE) #format into a matrix
# country2<-pred.labs.extract2[,1] 
# state2<-pred.labs.extract2[,2]
# #country.labs2<-c('AR', 'BR','CO' ,'EC', 'DR', 'GY', 'NC', 'PR') #MX, HR
# country.labs2<-countries
# cols.plot<- c('#e41a1c','#377eb8','#4daf4a','#984ea3', '#ff7f00','#a65628','#f781bf', 'grey')
# dim1<- max(1,floor(sqrt(length(countries))))
# dim2<- max(1,ceiling(sqrt(length(countries))))
# if(dim1+dim2<length(countries)){dim1=dim1+ 1}
# tiff(paste0(output_directory,'small multiples nobias.tiff'), width = 7, height = 4, units = "in",res=200)
# par(mfrow=c(dim1,dim2) , mar=c(2,2,0.6,0))
# for (i in 1:dim(preds.nobias.q.alt )[3]){
#   y.all<- t(exp(preds.nobias.q.alt[,,i]))
#   keep.obs<-as.data.frame(complete.cases(y.all))[,1]
#   y.all<-y.all[keep.obs,]
#   n.times<- nrow(y.all)
#   xx<- c(1:n.times, rev(1:n.times))
#   yy<- c(y.all[,1], rev(y.all[,3] ))
#   matplot(t(exp(preds.nobias.q.alt[,,i])), bty='l', type='l', col='white', lty=c(2,1,2), xlim=c(0,tot_time), ylim=c(0.5,2), axes=F)
#   polygon(xx, yy, col='gray90', border=NA)
#   points( y.all[,2], col=cols.plot[country2[i]], type='l')
#   abline(h=1, lty=2, col='black', lwd=1)
#   abline(h=c(0.8,1.2), lty=2, col='gray', lwd=0.5)
#   Axis(side=1, labels=FALSE, at=c(0,12,24,36,48,60), tck=-0.015,  col='gray')
#   Axis(side=2, labels=FALSE, at=c(0.5,1,2.0), tck=-0.005, col='gray')
#   if(state2[i]==1){
#     title(country.labs2[country2[i]], col.main=cols.plot[country2[i]] )
#     }
# }
# dev.off()
# 
# 
}

#################################
#################################
#COMBINE TOGETHER ESTIMATES FROM DIFFERENT AGE GROUPS
for(agegrp.select in c( '<2m','2-11m',  '12-23m', '24-59m')){
 print(agegrp.select)
   age_group <- agegrp.select # <2m, 2-11m, 2-23m, 2-59m, 12-23m, 24-59m
  output_directory<- paste0(dirname(getwd()), "/Results CP/",'A',"_",age_group, "_nat_MVN", '0','/')
  output_directory<-gsub("<2", "u2", output_directory)
  
  ds1<-readRDS( file=paste0(output_directory,"reg_mean_with_pooling CP nobias.rds"))
  lab1<-readRDS( file=paste0(output_directory, "country_labs.rds"))
  
  preds.q<- apply(ds1,c(2,3),  quantile, probs=c(0.025,0.5,0.975), na.rm=TRUE)
  dimnames(preds.q)[[3]]<-lab1
  agegrp.lab<-  gsub("<2", "u2", age_group)
  agegrp.lab<-  gsub("-", "_", agegrp.lab)
  
  assign(paste0("DF", agegrp.lab), preds.q)
}
#Peru

pr.u2<-t(DFu2m[,,'PAHO_pr'])
pr.2_11m<-t(DF2_11m[,,'PAHO_pr'])
pr.12_23m<-t(DF12_23m[,,'PAHO_pr'])
pr.24_59m<-t(DF24_59m[,,'PAHO_pr'])

tiff(paste0(output_directory,'pr.age.comp.tiff'), width=5, height=4 , units='in', res=200)
col.plot<-c('#e41a1c','#377eb8', '#4daf4a', '#984ea3')
matplot(((1:tot_time)-pre.vax.time),pr.u2, type='l', col=col.plot[1],bty='l', ylim=c(-1,0.5), lty=c(2,1,2), lwd=c(0.5,2,0.5))
matplot(((1:tot_time)-pre.vax.time),pr.2_11m, type='l', col=col.plot[2],bty='l', ylim=c(-1,0.5), lty=c(2,1,2), lwd=c(0.5,2,0.5),add=TRUE)
matplot(((1:tot_time)-pre.vax.time),pr.12_23m, type='l', col=col.plot[3],bty='l', ylim=c(-1,0.5), lty=c(2,1,2), lwd=c(0.5,2,0.5),add=TRUE)
matplot(((1:tot_time)-pre.vax.time),pr.24_59m, type='l', col=col.plot[4],bty='l', ylim=c(-1,0.5), lty=c(2,1,2), lwd=c(0.5,2,0.5),add=TRUE)
legend(-10,-0.4,legend=c('<2m','2-11m','12-23m','24-59m'), col=col.plot, lty=1,box.lty=0)
dev.off()