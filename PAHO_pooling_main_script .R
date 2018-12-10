
#Loop trough all of the age categories
#for(agegrp.select in c('<2m', '2-11m', '2-23m', '2-59m', '12-23m', '24-59m')){
rm(list=ls(all=TRUE))
for(agegrp.select in c( '<2m','2-11m', '2-23m', '2-59m', '12-23m', '24-59m')){
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
output_directory<- paste0(dirname(getwd()), "/Results/",hdi_level,"_",age_group, "_nat_MVN", max(subnational),'/')
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
countries<-countries[!nodata.country]
saveRDS(countries, file=paste0(output_directory, "country_labs.rds"))
