######## GROUP FOUR ELIZABETHKINGIA 


####dependencies
install.packages("triangle")
library(triangle)


########################### RISK ENDOTRACHEAL TUBE################################################

#####Concentration in tap water
conc_min<-1
conc_max<-5120
conc_avg<-12.25 

mL_on_hands<-0.5/4 ####per 1 hand fronts

#####transfer rate from hand to fomite accounts for 2 fronts of hands
transfer_low<-0.13
transfer_high<-0.38

#####area of endotracheal tube


#####area of fingertips
sa_fingertips_mean<-43.88 ###cm^2 one hand
sa_fingertips_sd<-1.707
hand_area_mean<- 365.3 ##one hand surface
hand_area_sd<-60.1
fraction_fingertips_hand<-rnorm(n,sa_fingertips_mean,sa_fingertips_sd)/rnorm(n,hand_area_mean,hand_area_sd)

###Area of tube
sa_tube_mean<-746.44 ### calculation(22.2+23)/2
sa_tube_sd<-1.860108 ####sqrt(1.1^2+1.5^2)

#####transfer rate from fomite to person
transfer_fomite<-1



####number of monte carlo samplings
n=10000



#######Common response variables
morbidity<-runif(n,0.19,0.78) ###need a range
mortality<-rtriangle(n,0,1,0)

####Dose calculation
dose<-rtriangle(n,conc_min,conc_max,conc_avg)/100*mL_on_hands*fraction_fingertips_hand*runif(n,transfer_low,transfer_high)*transfer_fomite 


#placeholder<-c(0.32,1,0.18,0,.33,0.5,0.29,0.65,0,0.2,0,.25,0,.23,0.31)
#hist(placeholder,breaks=20)
#hist(morbidity)


#######RESPONSE BETA POISSON (contact lense)
alpha<-1.9E-01	
N50<-1.85E+04
Pinf_bp<-1-((1+dose*((2^(1/alpha)-1)/N50)))^-alpha
Pill_bp<-Pinf_bp*morbidity
Pdeath_bp<-Pill_bp*mortality

########RESPONSE EXPONENTIAL (injection)

##first need to calculate sd from percentile data provided in the wiki
#p95<-log(1.48E-04)#from wiki

mean_k <- 0.000105 ##from wiki
top_k <- 0.000157
z_crit <- 1.95996398454005
ln_k <- log(mean_k)
#ln_k #this verifies mean k
ln_topk<- log(top_k)   
sd_k <- (ln_topk-ln_k)/z_crit
#sd_k #this verifies sd for k

##k is a log normal distribution using parameters calculated above
k<-rlnorm(n,ln_k,sd_k)
Pinf_e<-1-exp(-k*dose)
Pill_e<-Pinf_e*morbidity
Pdeath_e<-Pill_e*mortality

number_intubated_year<-1701440


#######Visualizations
endo_risk<-data.frame(Pinf_bp,Pill_bp,Pdeath_bp,Pinf_e,Pill_e,Pdeath_e)
#####Per year
Pinf_bp_yr<-number_intubated_year*Pinf_bp
Pill_bp_yr<-number_intubated_year*Pill_bp
Pdeath_bp_yr<-number_intubated_year*Pdeath_bp
Pinf_e_yr<-number_intubated_year*Pinf_e
Pill_e_yr<-number_intubated_year*Pill_e
Pdeath_e_yr<-number_intubated_year*Pdeath_e

endo_risk_cases_year<-data.frame(Pinf_bp_yr,Pill_bp_yr,Pdeath_bp_yr,Pinf_e_yr,Pill_e_yr,Pdeath_e_yr)
#####images for export
####Risk images

png("Risk_Endotracheal_tube_Beta_Poisson_Exponential.png",width = 8, height = 8, units = 'in', res = 1200)  
par(mfrow=c(1,1))
boxplot(endo_risk, main=
          "Risk of infection, illness, and death from E. anophelis-contaminated
           endotracheal tubes in immunocompromised individuals separated 
           by dose response model", log="y",
        ylab="Probability",
        xlab="Beta Poisson Model         |         Exponential Model",
        names= c("Infection","Illness","Mortality","Infection","Illness","Mortality"),
        col = c("turquoise2","turquoise3","turquoise4","violetred2","violetred3","violetred4"))
abline(v=3.5)
dev.off()


### PER year
png("Cases_Endotracheal_tube_Beta_Poisson_Exponential.png",width = 8, height = 8, units = 'in', res = 1200)  
par(mfrow=c(1,1))
boxplot(endo_risk_cases_year,
        cex.main=0.985,main=
          "Cases of infection, illness, and death caused by E. anophelis on endotracheal tubes
            in immunocompromised individuals per year in the United States separated 
               by dose response model",
        ylab="Potential Cases/Year",
        #ylim=c(0,150),
        xlab="Beta Poisson Model         |         Exponential Model",
        names= c("Infection","Illness","Mortality","Infection","Illness","Mortality"),
        col = c("turquoise2","turquoise3","turquoise4","violetred2","violetred3","violetred4"))
abline(v=3.5)
dev.off()
###mortality only

png("Mortality_Endotracheal_tube_Beta_Poisson_Exponential.png",width = 8, height = 8, units = 'in', res = 1200)  
par(mfrow=c(1,1))
boxplot(endo_risk_cases_year$Pdeath_bp_yr,endo_risk_cases_year$Pdeath_e_yr,main=
          "Mortality cases resulting from E. anophelis contaminated endotracheal 
          tubes immunocompromised individuals per year in the United States
          separated by dose response model*",
         ylab="Potential Cases/Year",
          yaxt = "none",
        names= c("Mortality (Beta Poisson Model)","Mortality (Exponential Model)"),
        col = c("turquoise4","violetred4"),
          ylim=c(0,30),
        outline=FALSE)
axis(2, seq(0,30,1),cex.axis=1,las=2, font=1)
mtext(side=1, line=2, "***Outliers beyond whisker go beyond graph bounds but are excluded for visualization", font=1,cex=0.8)
dev.off()

####Dose images

png("Dose_Endotracheal_tube.png",width = 8, height = 8, units = 'in', res = 1200)  
par(mfrow=c(1,1))
hist(dose,xlab="Dose in CFUs",main=
       "Dose of E. anophelis in an contaminated endotracheal tube scenario",
        col="skyblue1",
        xlim=c(0,0.6),
        breaks=50)
rug(dose)
dev.off()



############################ RISK FEMALE CATHETER ################################################

############################ RISK MALE CATHETER ################################################


