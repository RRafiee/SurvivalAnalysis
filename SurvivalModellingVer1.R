##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Survival modelling script written by Dr Reza Rafiee
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script gets a csv file including all variables for doing survival modelling analysis

install_all_packages_automatic <- function(x) {
  x <- as.character(substitute(x))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    #update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

# Install all required packages as well as dependencies
install_all_packages_automatic(pec)
install_all_packages_automatic(pROC)
install_all_packages_automatic(survival)
install_all_packages_automatic(survivalROC)

# Loading all survival and dependent libraries
library(survival)
library(survivalROC)
library(pec)
library(pROC)


Pathfolder <- "~/survivalmodelling/"  # initialise the path of the csv file
csvfilename <- "survivaldata.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)

# Checking the prorpotionality of hazard
var1_ph <- cox.zph(coxph(Surv(OS_Time, OS_YN) ~ Variable1, data = datafileObject))
var1_ph_pvalue <- var1_ph$table[3] # p-value
#            rho    chisq    p
#Variable1  0.281  1.56    0.212   # p-value: 0.212, is not significant ==> doesn't reject the null hypothesis (Pass the test for the proportionality of hazard)
var2_ph <- cox.zph(coxph(Surv(OS_Time, OS_YN) ~ Variable2, data = datafileObject))
var2_ph_pvalue <- var2_ph$table[3] # p-value
#            rho   chisq     p
#Variable2 -0.0872 0.134   0.714   # p-value: 0.714, is not significant ==> doesn't reject the null hypothesis (Pass the test for the proportionality of hazard)

par(mfrow=c(2,2))

plot(main= "proportional hazard test", cox.zph(coxph(Surv(OS_Time, OS_YN) ~ Variable1, data = datafileObject))) # variable1
legend(bty="n","topleft",lty=0.1, col="blue", paste("p-value: ",round(var1_ph_pvalue,digits = 3)))
plot(main= "proportional hazard test", cox.zph(coxph(Surv(OS_Time, OS_YN) ~ Variable2, data = datafileObject))) # variable2
legend(bty="n","topleft",lty=0.1, col="blue", paste("p-value: ",round(var2_ph_pvalue,digits = 3)))

# Building a survival model
Survival_model <- coxph(Surv(OS_Time, OS_YN) ~ Variable1 +  Variable2, data = datafileObject)
summary(Survival_model)

# summary of the Survival_model:
# n= 57, number of events= 19 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# Variable1 1.6213    5.0599   0.4777 3.394 0.000689 ***
# Variable2 1.7622    5.8252   0.4734 3.722 0.000197 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Variable1     5.060     0.1976     1.984     12.91
# Variable2     5.825     0.1717     2.303     14.73
# 
# Concordance= 0.814  (se = 0.058 )
# Rsquare= 0.361   (max possible= 0.917 )
# Likelihood ratio test= 25.51  on 2 df,   p=2.892e-06
# Wald test            = 25.81  on 2 df,   p=2.489e-06
# Score (logrank) test = 34.26  on 2 df,   p=3.636e-08
# Both variables are significant
## Predict survival at 5 years (for the input csv file)
survival_probability <- predictSurvProb(Survival_model, newdata=datafileObject, times = 5)
## Examine predicted 5 year survival for this csv file 
boxplot(survival_probability, outpch=NA, main="Survival ROC at 5 years", ylab="csv cohort 5 year survival probability", ylim=c(0,1))
stripchart(method="jitter", add=TRUE, vertical=TRUE, survival_probability, col="blue", pch=19)
abline(h=c(0.6), col="red")


AUC_survivalmodel <- survivalROC(Stime = datafileObject$OS_Time,  
                        status = datafileObject$OS_YN,      
                        marker = 1-survival_probability,  
                        predict.time =  5, method="KM") 

AUC_survivalmodel$AUC 
# 0.8543472


plot(main="Survival ROC at 5 years", AUC_survivalmodel$FP, AUC_survivalmodel$TP, type="l", col="blue",xlim=c(0,1), ylim=c(0,1),   
     xlab="FP", ylab="TP", lwd=2, axes=F)
axis(1)
axis(2, las=2)
abline(col="grey", c(0,1))


## Simplified the risk model
simplifiedmodel <- ifelse(survival_probability < 0.6, 1,0)
AUC_simplifiedmodel <- survivalROC(Stime = datafileObject$OS_Time,  
                             status = datafileObject$OS_YN,      
                             marker = simplifiedmodel,    
                             predict.time =  5, method="KM") 

AUC_simplifiedmodel$AUC # 0.8374794 

## Predict survival at 5 years
lines(col="red",AUC_simplifiedmodel$FP, AUC_simplifiedmodel$TP, type="l" )


legend(bty="n", "bottomright",lty=1, col=c("blue","red"), c(paste("Cox Model - ", round(AUC_survivalmodel$AUC, digits = 3),sep = ""),
                                                            paste("Simplified Model - ", round(AUC_simplifiedmodel$AUC,digits = 3),sep = "")))

temp_probability_survival <- data.frame(survival_probability, simplifiedmodel, datafileObject)
temp_probability_survival[order(temp_probability_survival$survival_probability),]
write.table(temp_probability_survival[order(temp_probability_survival$survival_probability),], file="Survival_probabilities.csv", sep=",")


#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
