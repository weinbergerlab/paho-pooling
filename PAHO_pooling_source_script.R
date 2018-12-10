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
  
  index.post<-which(time>=(eval_period[1] %m-% months(pre.vax.time)) & time<=eval_period[2])  
  if (length(index.post)>tot_time){ index.post<-index.post[1:tot_time] } 
  
  #Test whether RR differs by month
  month<-as.factor(month(time))
  post.index.p12<- time>=(eval_period[1] %m+% months(12))
  log.rr.median<-log_rr_q[,2,, drop=FALSE]
  log.rr.median<-matrix(log.rr.median[,1,],nrow=nrow(log.rr.median)) 
  log.rr.median.spl<-split(t(log.rr.median),1:ncol(log.rr.median)) 
  seas.func1<-function(ds){
   mod1<-lm(ds~month*post.index.p12)
   summary(mod1)
   return(summary(mod1))
  }
  mods<-lapply(log.rr.median.spl, seas.func1 )
  names(mods)<-dimnames(log_rr_q)[[3]]
  print(mods)
 # agegrplab<-agegrp.select
  #if(agegrp.select=='<2m'){agegrplab<-'u2m'}
  saveRDS(mods, file=paste0(output_directory,'_', country,"Seasmods.rds"))
}

  
  