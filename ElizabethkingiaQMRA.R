######## GROUP FOUR ELIZABETHKINGIA 


####dependencies
install.packages("triangle")
library(triangle)


#### DOSE

n=10000
dose<-rtriangle(n,1,5120,12.25)/100*0.5*runif(n,0.065,0.19)*1 
hist(dose)
Pinf