######## GROUP FOUR ELIZABETHKINGIA 


####dependencies
install.packages("triangle")
library(triangle)


########################### RISK ENDOTRACHEAL TUBE################################################

#####Concentration in tap water
conc_min<-1
conc_max<-5120
conc_avg<-12.25 

mL_on_hands<-0.5 ####per pair of hands

#####transfer rate from hand to fomite accounts for half of hand surfaces
transfer_low<-0.065
transfer_high<-0.19


#####area of endotracheal tube


#####area of fingertips


#####transfer rate from fomite to person
transfer_fomite<-1


####number of monte carlo samplings
n=10000

dose<-rtriangle(n,conc_min,conc_max,conc_avg)/100*mL_on_hands*runif(n,transfer_low,transfer_high)*transfer_fomite 

#######Common response variables
morbidity<-0.67 ###need a range
mortality<-runif(n,0.23,0.52)

#######RESPONSE BETA POISSON (contact lense)
alpha<-1.9E-01
N50<-1.85E+04
Pinf_bp<-1-((1+dose*((2^(1/alpha)-1)/N50)))^-alpha
Pill_bp<-Pinf_bp*morbidity
Pdeath_bp<-Pill_bp*mortality

########RESPONSE EXPONENTIAL (injection)

##first need to calculate sd from percentile data provided in the wiki
p95<-log(1.48E-04)
meank<-log(1.05E-04)
zp95<-1.645
sd<-(p95-meank)/zp95
##k is a log normal distribution using parameters calculated above
k<-rlnorm(10000,meank,sd)
Pinf_e<-1-exp(-k*dose)
Pill_e<-Pinf_e*morbidity
Pdeath_e<-Pill_e*mortality


#######Visualizations
endo_risk<-data.frame(Pinf_bp,Pill_bp,Pdeath_bp,Pinf_e,Pill_e,Pdeath_e)

png("Risk_Endotracheal_tube_Beta_Poisson_Exponential.png")  
par(mfrow=c(1,1))
boxplot(endo_risk, main=
          "Risk of infection, illness, and death separated 
        by dose response model for endotracheal tubes in
        immunocompromised individuals", log="y")
dev.off()
boxplot(endo_risk, main=
          "Risk of infection, illness, and death separated 
        by dose response model for endotracheal tubes in
        immunocompromised individuals", log="y")

hist(dose,xlab="dose in CFUs",main="Dose of E. anophelis in an endotracheal tube scenario")

png("Dose_Endotracheal_tube.png")  
par(mfrow=c(1,1))
hist(dose,xlab="dose in CFUs",main="Dose of E. anophelis in an endotracheal tube scenario")
dev.off()



############################ RISK FEMALE CATHETER ################################################

############################ RISK MALE CATHETER ################################################


