N.countries=length(countries)
log.rr.compile <-vector("list", length(countries)) 
log.rr.prec.mat.all<-vector("list", length(countries)) 
N.states<-rep(NA, N.countries)

for (c in 1:N.countries){
  country<-countries[c]
  print(country)
  
  
  #####################################################################
  #SET COUNTRY_SPECIFIC PARAMETERS--set eval period as intro date-1 month
  if (grepl("PAHO_ar", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2012-01-01'), as.Date('2015-12-01'))
    keep.grp=paste0('ar ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_br", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2010-03-01'), as.Date('2015-12-01'))
    keep.grp=paste0('br ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_co", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2011-11-01'), as.Date('2015-12-01'))
    keep.grp=paste0('co ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_dr", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2013-09-01'), as.Date('2015-12-01'))
    keep.grp=paste0('dr ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_ec", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2010-08-01'), as.Date('2015-12-01'))
    keep.grp=paste0('ec ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_gy", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2011-01-01'), as.Date('2013-12-01'))
    keep.grp=paste0('gy ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_hr", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2011-01-01'), as.Date('2015-12-01'))
    keep.grp=paste0('hr ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_mx", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2008-01-01'), as.Date('2015-12-01'))
    keep.grp=paste0('mx ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_nc", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2012-01-01'), as.Date('2015-12-01'))
    keep.grp=paste0('nc ', age_group, ' ', hdi_level)
  } else if (grepl("PAHO_pr", country, fixed=TRUE)) {
    eval_period <- c(as.Date('2009-08-01'), as.Date('2014-12-01'))
    keep.grp=paste0('pr ', age_group, ' ', hdi_level)
  } 
  ####################################################################
  
  #Import the data for each country and extract the relevant sub groups and time periods
  input_directory <- paste(dirname(getwd()),'Data',country, sep='/')
  log_rr_q<-readRDS(file=paste0(input_directory,'/', country, "_log_rr_quantiles_best.rds"))
  log_rr_prec<-readRDS(file=paste0(input_directory,'/', country, "_log_rr_best_t_samples.prec.rds"))
  time<-as.Date(as.numeric(dimnames(log_rr_q)[[1]]), origin="1970-01-01" )  
  
  index.post<-which(time>=eval_period[1] & time<=eval_period[2])  
  if (length(index.post)>max.time.points){ index.post<-index.post[1:max.time.points] } 
  
  if(subnational[c]==0){  #National-level only
    # limit by age group
    keep.index<-which(dimnames(log_rr_q)[[3]]== (keep.grp))
    log_rr_q<-log_rr_q[index.post,,keep.index]
    log_rr_prec<-log_rr_prec[index.post,index.post,keep.index]
    
    log_rr_q<-array(log_rr_q,dim=c(nrow(log_rr_q),ncol(log_rr_q),1))
    log_rr_prec<-array(log_rr_prec,dim=c(nrow(log_rr_prec),ncol(log_rr_prec),1))
  } else {  
    select.obs<-grep(paste0("^",keep.grp), dimnames(log_rr_q)[[3]])
    exclude.obs2<-which(apply(log_rr_prec,3,function(x)  sum(is.nan(x)))>0)  #Which states have no NaN in covariance matrix
    
    select.obs.state<-setdiff(select.obs,exclude.obs2)      
    log_rr_q <-log_rr_q[index.post,,select.obs.state]
    #log_rr_sd <-log_rr_sd[,index.post,]
    log_rr_prec<-log_rr_prec[index.post,index.post,select.obs.state]
  }
  N.states[c]<-dim(log_rr_q)[3]
  log.rr.compile[[c]]<-log_rr_q
  log.rr.prec.mat.all[[c]]<-log_rr_prec
}

##Combine together estimates from each country into a single array
log_rr_q_all<-array(NA, dim=c(max.time.points,max(N.states),N.countries))
state.labels<-array(NA, dim=dim(log_rr_q_all)[c(2:3)]   )
log_rr_prec_all<-array(NA, dim=c(max.time.points,max.time.points,max(N.states),N.countries))
ts.length<-rep(NA, times=N.countries)
for(i in 1:N.countries){
  for(j in 1:N.states[i]){
    log_rr_q_all[1:nrow(log.rr.compile[[i]]),j,i]<-log.rr.compile[[i]][,2,j] #Extract the median
    log_rr_prec_all[1:dim(log.rr.prec.mat.all[[i]])[1],1:dim(log.rr.prec.mat.all[[i]])[1],j,i]<-log.rr.prec.mat.all[[i]][,,j]
  }
  ts.length[i]<-nrow(log.rr.compile[[i]])
  if(N.states[i]>1){
    state.labels[1:N.states[i] ,i]<-dimnames(log.rr.compile[[i]])[[3]]
  }else{
    state.labels[1,i]<-'national'
  }
}
N.time.series <- sum(N.states)

#PREPARE SPLINE BASIS using gam()
set.seed(15785)
N.knots <- 4  #N-1 internal knots
log.rr.q.extr<-log_rr_q_all[,,1]
sp.df<-as.data.frame(log.rr.q.extr)
names(sp.df)<-paste0('state', 1:ncol(sp.df))
sp.df$time=1:nrow(sp.df)
sp.m<-melt(sp.df ,id='time'  )
sp.c<- acast(sp.m, variable~time)
t.index<-1:max.time.points
out<-rnorm(n=max.time.points)
ds.sp<-cbind.data.frame(out,t.index)

spl.t.std<-bs(1:max.time.points, degree=3,df=N.knots, intercept=FALSE) 
rowSums(spl.t.std)
matplot(spl.t.std, type='l')
spl.t.std<-cbind(rep(1, nrow(spl.t.std)), spl.t.std)

ts.length_mat<-matrix(NA, nrow=N.countries, ncol=max(N.states))
for(i in 1:N.countries){
  ts.length_mat[i,1:N.states[i]]<-ts.length[i]
}
identity<-diag(nrow(spl.t.std))
ts.length.vec<-as.vector(t(ts.length_mat))[!is.na(as.vector(t(ts.length_mat)))]

#Matrices to input to JAGS
#I_Sigma<-replicate( N.countries, diag(p) )
I_Sigma<-diag(p)  #For p=N vcovariates of intercept and slope
I_Omega<-diag(q)  # q= number of predictors of the p slopes 
m<-1  #m=number of country-level covariates
w<-matrix(NA, nrow=N.countries, ncol=m) 
for(i in 1:N.countries){
  for(k in 1:m){
    w[i,k]<-1
  }
}
