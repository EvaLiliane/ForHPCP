rm(list=ls(all=T))
library(zoo)
library(longitudinal)
library(foreign)
library(nlme) 
library(plyr) 
library(lattice)
#library(lme4)

library(mgcv) 
library(ggplot2) 
library(reshape2)
setwd("/home/evaliliane/Documents/PhD/Codes")

# Read in the two children dataset
#addata <- read.csv("/home/evaliliane/Documents/PhD/Codes/NewData/CD4Cat_Adults2015-04-20.csv")
chdata <- read.csv("/home/evaliliane/Documents/PhD/Codes/NewData/CD4Cat_Children2015-04-20.csv")
# str(chdata)
# rel.columns <- c("patient", "diff", "zscore", "lab_v","height","weight","viral", "fhv_who_stage", "gender")
# ddat2red <- chdata[chdata$base == 1, c("patient", "height","weight","viral", "fhv_stage_who", "gender")]
# names(ddat2red) <- c("patient", "baseheight","baseweight","baseviral", "basefhv_stage_who", "basegender")
                   
                   

# Replacing small neg zscores by -12
chdata$zscore[chdata$zscore < -10] <- -10

# Removing row with TOFU > 15
chdata <- chdata[chdata$diff < 5476,]
dat <- chdata[,c(2,35,44,48,50)]
dat <- dat[!(is.na(dat$z.categ)),]

# Select a subgroup of people
useSubset <- T
subselect <- function(addata,n){ 
  num.to.do <- n ## to save computing time when playing with analyses
  pids.to.run1 <- sample(unique(addata$patient), num.to.do)
  addata <- addata[addata$patient %in% pids.to.run1,] 
  return(addata)
}
datt <- subselect(dat, 1000)

#subchdata <- dat


currentDate <- Sys.Date()

# ##############################################################################################
# ################      Model Building         #################################################
# ##############################################################################################


# ===================== Children per CD4 Categories ==========================

testchdata <- datt #[addata$suppress == 1,]      #subselect(addata,2000)
nn <- nlevels(testchdata$z.categ)
mydataa <- groupedData(lab_v ~ diff | z.categ/patient, data = testchdata,order.groups=F) # Only On-ART  period

# Order the data and get intial values for the parameter estimates
mydataa <- mydataa[order(mydataa$z.categ, mydataa$patient,mydataa$diff),]  
# model1<-nls(zscore ~ SSasymp(diff,Asym,R0,lrc), data=mydataa,na.action="na.omit",
#               start=c(Asym=-2,R0=-3,lrc=-5), control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/4096,
#                                                                    printEval = FALSE, warnOnly = FALSE)) 



#Individual fits for Old Model
model1.lis <- nlsList(zscore ~ SSasymp(diff,Asym,R0,lrc)|z.categ, data=mydataa,
                start=c(Asym=-2,R0=-4,lrc=-1),
                pool = FALSE)

model2.lis <- nlsList(zscore ~ SSasymp(diff,Asym,R0,lrc)|z.categ, data=mydataa,
                      start=c(Asym=-2,R0=-4,lrc=-1))

model3.lis <- nlsList(zscore ~ SSasymp(diff,Asym,R0,lrc)|patient, data=mydataa,
                      start=c(Asym=-2,R0=-4,lrc=-1))
model4.lis <- nlsList(zscore ~ SSasymp(diff,Asym,R0,lrc)|patient, data=mydataa,
                      start=c(Asym=-1,R0=-3,lrc=-1))
model5.lis <- nlsList(zscore ~ SSasymp(diff,Asym,R0,lrc)|patient, data=mydataa,
                      start=c(Asym=-1.8,R0=-3.9,lrc=-1.7))
# Comments : Pool does not help

# Plot confidence intervals
plot(intervals(modeltest), layout = c(5,1), na.rm = TRUE)
plot(intervals(model5.lis), layout = c(3,1), na.rm = TRUE)
plot(intervals(model4.lis), layout = c(3,1), na.rm = TRUE)
plot(intervals(model3.lis), layout = c(3,1))
plot(intervals(model2.lis), layout = c(3,1))
#anova(model2.lis,model1.lis)
# Results show that random effects are needed for all three parameters.
cor(coef(model3.lis), use = "pairwise.complete.obs")
sds <- sappply(model2.lis, sigmaHat)

tt <- coef(model5.lis)
tt <- cbind(patient = rownames(tt), tt)
convt <- tt$patient[!(is.na(tt$R0))]

pp <- coef(model4.lis)
pp <- cbind(patient = rownames(pp), pp)
convp <- pp$patient[!(is.na(pp$R0))]

convdata <- mydataa[mydataa$patient %in% convt,]
myconv <- groupedData(zscore ~ diff | z.categ/patient, data = convdata,order.groups=F)
plot(myconv)
notconvdata <- mydataa[!(mydataa$patient %in% convp),]
mynotconv <- groupedData(zscore ~ diff | z.categ/patient, data = notconvdata,order.groups=F)
plot(mynotconv[1:100,])


# Look at the residual per z.categ
plot(model2.lis, z.categ ~ resid(.), abline = 0)

# Make individual plots
plot(mydataa[550:600,])
plot(mydataa, outer = ~ diff * z.categ)
plot(augPred(model2.lis[580:600,]))
# Individual sum of LogLik
logLik(model2.lis)  # Does not work in presence of NULL fits

# Scatter plot of the estimated parameters relationships
pairs(model2.lis)
pairs(model3.lis)


cof <- coef(summary(model1))
init = cof[1:3] #c(Asym = 450, R0 = 380, lrc = -6)

# # ======================================================================
# # Single level models
# #Random effect defined for CD4 categories level
#model2a <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
#               data = mydataa, 
#               na.action="na.omit",
#               random = Asym + R0 ~ 1|cd4a.categ, 
#               fixed = list(Asym ~ 1, R0 ~ 1, lrc ~ 1),                        
#               start= init, verbose = FALSE )
#
##Random effect defined for patient level
#model2b <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
#               data = mydataa, 
#               na.action="na.omit",
#               random = Asym + R0 ~ 1|patient, 
#               fixed = list(Asym ~ 1, R0 ~ 1, lrc ~ 1),                        
#               start= init, verbose = FALSE )
#
# # Look at both outputs.
# summary(model2a)
# #          StdDev    Corr 
# # Asym     0.5431035 Asym 
# # R0       2.4207631 0.989
# # Residual 1.9394735 
# summary(model2b)
# #          StdDev   Corr 
# # Asym     1.748278 Asym 
# # R0       2.847446 0.482
# # Residual 1.369344

# ======================================================================
# Multiple levels models
# Random effect defined for CD4 categories and patient levels
model3 <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
               data = mydataa, 
               na.action="na.omit",
               random = Asym + R0 ~ 1|z.categ/patient,
               fixed = list(Asym ~ 1, R0 ~ 1, lrc ~ 1),                        
               start= init, verbose = FALSE ) 
print("Model 3 with  fixed effect and a two levels random effect.")
summary(model3)
plot(model3)
# High correlation between the parametres Asym & R0
# Plot shows an increase in the variability of the residuals

# modeltest <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
#                data = mydataa, 
#                na.action="na.omit",
#                random = Asym + R0 + lrc ~ 1 , #|z.categ/patient, lrc ~ 1 ),#|z.categ/patient),
#                fixed = list(Asym ~ 1, R0 ~ 1, lrc ~ 1),                        
#                start= init, verbose = TRUE ) 
# 
# Error in nlme.formula(zscore ~ SSasymp(diff, Asym, R0, lrc), data = mydataa,  : 
#                         Singularity in backsolve at level 0, block 1
#                       In addition: Warning message:
#                         In nlme.formula(zscore ~ SSasymp(diff, Asym, R0, lrc), data = mydataa,  :
#                                           Singular precision matrix in level -1, block 4
# --------------------------------------
# Model 3 corrected for independant random effects on both levels
model4 <- update(model3, random = list(z.categ = pdDiag(Asym + R0 ~ 1),patient = pdDiag(Asym + R0 ~ 1)))
print("Model3 corrected for independant random effects on both levels.")
summary(model4)
#Compare bothe model 3 & model4
anova(model3, model4)
# Results: According to the AIC values, there is no improvement from model3 to model4. However, addressing the
#          correlation problem between the parameters seem of higher priority than reducing the AIC 

# ----------------------------
# Inclusion of a covarience matrix in Model 4 
# !!!!!!! NOT WORKING
# Error in Initialize.corARMA(X[[2L]], ...) : 
#   covariate must have unique values within groups for "corARMA" objects
# model5  <- update(model4, corr = corARMA(c(0.5,0.5),form = ~ diff |z.categ/patient , p=1,q=1))
# # print("Inclusion of a covarience matrix in Model 4.")
# summary(model5)
# # Compare model 4 & model 5
# anova(model4, model5)
# results: 


# ----------------------------
# Model 5 corrected for heteroscedasticity
model66 <- update(model3, weights = varPower())
model6 <- update(model3, weights = varIdent( 0.2, ~ 1|z.categ))
print("Model5 corrected for heteroscedasticity.")
summary(model6)
# Compare model 5 & model 6
anova(model4, model6)
# results: 


# -------------------------------------------------------
# Try Univariate models for different covariates
# Age at baseline
new.cof <- model3$coeff
new.init1 <- c(as.numeric(new.cof$fixed[1:3]),0,0,0)
model6age1 <- update(model3, fixed = list(Asym ~ age, R0 ~ age, lrc ~ age), start = new.init1) 


new.init <- c(as.numeric(new.cof$fixed[1:3]),0,0)
model6age2 <- update(model3, fixed = list(Asym ~ age, R0 ~ age, lrc ~ 1), start = new.init)
print("Model6 with Age covariate.")

# Compare model 6 & model.Age
#anova(model3, model6age)
anova(model6age1, model6age2)
summary(model6age1)
# Results: There is a significant improvement when including Age as a confounder for the Fixed parameters

# # Gender
model6gender <- update(model3, fixed = list(Asym + R0 ~ gender, lrc ~ 1), start = new.init)
print("Model6 with Gender covariate.")
summary(model6gender)
# Compare model 6 & model.gender
anova(model3, model6gender)
# # Results: 
# 
# # Baseline WHO stage
# model6who <- update(model3, fixed = list(Asym + R0 ~ fhv_stage_who, lrc ~ 1), start = new.init)
# print("Model6 with Who stage covariate.")
# summary(model6who)
# # Compare model 6 & model.who
# anova(model3, model6who) # Not comparable bcz of missing values in WHO stage variable
# Results:

# # Baseline BMI 
# model6bmi <- update(model3, fixed = list(Asym + R0 ~ bmi, lrc ~ 1), start = new.init)
# print("Model6 with bmi covariate.")
# summary(model6bmi)
# # Compare model 6 & model.bmi
# anova(model3, model6bmi)
# # Results:

# # Baseline RNA 
# model6rna <- update(model3, fixed = list(Asym + R0 ~ viral, lrc ~ 1), start = new.init)
# print("Model6 with RNA covariate.")
# summary(model6rna)
# # Compare model 6 & model.rna
# anova(model3, model6rna)
# # Results:

# # Baseline OIs 
# model6rna <- update(model3, fixed = list(Asym + R0 ~ oi, lrc ~ 1), start = new.init)
# print("Model6 with Ois covariate.")
# summary(model6oi)
# # Compare model 6 & model.oi
# anova(model3, model6oi)
# # Results:

# -------------------------------------------------------
# Try full model of covariates
new.init <- c(as.numeric(new.cof$fixed[1:3]),0,0,0,0,0)
model6two <- update(model3, fixed = list(Asym ~ age + gender, R0 ~ age + gender, lrc ~ age), start = new.init) 
print("Model6 with Age and gender covariates.")
summary(model6two)
# Compare model one cov & model with two 
anova(model6age, model6two)
# results: No improvement by correcting for gender

# EXIT

# 
#  model33 <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
#                data = mydataa, #we want this to be mydataa right?  not mydata?
#                na.action="na.omit",
#                random = Asym + R0 ~ 1|cd4a.categ,
#                fixed = list(Asym ~ age + gender, R0 ~ age + gender, lrc ~ 1),#,+ fhv_stage_who + gender + suppress + weight + height, R0 ~ 1, lrc ~ 1),                        
#                start= c(Asym=-2,5,1,R0=-4,5,1,lrc=-5), verbose = FALSE )
#  
# # ================================ Asy .vs. Int ===============================
# r <- ranef(model3)
# f <- fixef(model3)
# r<- c( -2.840674, -1.049151, -2.377164, -1.908371, -1.602952)
# i <- c( -8.298243, -0.790337, -4.767319, -3.075800, -2.160361)
# 
# png("Output/Asym_Int.png", w=480,h = 480)
# plot(i,r, xlab = "Baseline z-score",pch = 15, ylab = "Long term z-score", col=1:6, main = "Relation between Asy & Int")
# #legend("toplef", levels(testchdata)[-6], col = 1:6, lty = 1, cex=0.7)
# dev.off()
# 
# # =============================================================================
# 
# barplot(r$Asym, names.arg = rownames(r))
# barplot(r$R0, names.arg = rownames(r))
# png("Output/Rebound_Int-Asym.png", w=480,h = 480)
# barplot(r$Asym - r$R0, names.arg = rownames(r))
# dev.off
# # =============================================================================
# 
# # Plot residualss
# png("Output/Residuals.png", w=480,h = 480)
# plot(fitted(model33),residuals(model33),main="Residuals vs Fitted", cex=0.5, xlab = "Predictions", ylab= "Standardized Residuals")
# dev.off()
# 
# # ##############################################################################
# 
# 
# # ===================== Children per CD4 Categories ==========================
# 
# testchdata <- modchdata #[addata$suppress == 1,]      #subselect(addata,2000)
# # #testchdata <- groupedData(zscore ~ diff | patient, data = testchdata,order.groups=F)
# # 
# # 
# nn <- nlevels(testchdata$cd4a.categ)
# mydataa <- groupedData(zscore ~ diff |patient, data = testchdata,order.groups=F) # Only On-ART  period
# mydataa <- mydataa[order(mydataa$patient,mydataa$diff),]  
# #print(paste0(length(unique(mydataa$patient)),'  ',sd(mydataa$lab_v), " patients count, std of CD4 counts, used for LME models"))
# model1<-nls(zscore ~ SSasymp(diff,Asym,R0,lrc), data=mydataa,na.action="na.omit",
#             start=c(Asym=-2,R0=-4,lrc=-5), control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/4096,
#                                                                  printEval = FALSE, warnOnly = FALSE))
# cof <- coef(summary(model1))
# #print(paste("Starting with category : ", curr.cat))
# print(cof)
# init = cof[1:3] #c(Asym = 450, R0 = 380, lrc = -6)
# 
# model1 <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
#                data = mydataa, #we want this to be mydataa right?  not mydata?
#                na.action="na.omit",
#                random = Asym + R0 ~ 1,
#                fixed = list(Asym ~ 1, R0 ~ 1, lrc ~ 1),#,+ fhv_stage_who + gender + suppress + weight + height, R0 ~ 1, lrc ~ 1),                        
#                start= c(Asym=-2,R0=-4,lrc=-5), verbose = FALSE ) 
# 
# 
# model2 <- nlme(zscore ~ SSasymp(diff,Asym,R0,lrc),
#                data = mydataa, #we want this to be mydataa right?  not mydata?
#                na.action="na.omit",
#                random = Asym + R0 ~ 1,
#                fixed = list(Asym ~ age + gender, R0 ~ age + gender, lrc ~ 1),#,+ fhv_stage_who + gender + suppress + weight + height, R0 ~ 1, lrc ~ 1),                        
#                start= c(Asym=-2,5,1,R0=-4,5,1,lrc=-5), verbose = FALSE ) 
# 
# # ================================ Asy .vs. Int ===============================
# r <- ranef(model2)
# f <- fixef(model2)
# plot(r$R0,r$Asym, xlab = "Baseline z-score", ylab = "Long term z-score", main = "Relation between Asy & Int")
# legend("toplef", rownames(r), lty = 1, cex=0.7)
# 
# # =============================================================================
# 
# barplot(r$Asym, names.arg = rownames(r))
# barplot(r$R0, names.arg = rownames(r))
# barplot(r$Asym - r$R0, names.arg = rownames(r))
# 
# # =============================================================================
# #plot(model33, resid(., type = "p") ~ fitted(.) | cd4a.categ, abline = 0)
# plot(fitted(model33),residuals(model33),main="Residuals vs Fitted", cex=0.5)
# plot(model33) # Same as above. This looks better
# 
# ####################################################
# # COMMENTS:
# # By comparing all nlme models, model5 sems to be the best fit to our data.
# # For the model formula, see refs Beaudrap (2008) and Lewis et al. (2011) in our Mendeley folder
# # KR: this formula is the same as the default SSasym. The only difference between model 4 and 5 is that
# # model 5 does not have a random R0 (intercept) which seems like it would be important in our
# # model since people start ART at many different starting CD4 counts and this has been shown to be
# # fairly important in determining what CD4 count people level off at
# 
# 
# ## Plot Model results - working
# modelpred.c1 <-as.numeric(predict(model2))
# modelpred.c2 <-labels(predict(model2))
# modelpred.c3 <-mydataa$diff
# modelpred <- data.frame(modelpred.c1,modelpred.c2, modelpred.c3)
# 
# nam <- c("predictions","Ids","years")
# names(modelpred) <- nam
# 
# # ====================================================== Functions
# 
# # Defined functions
# med_CI <- function (x, ci = 0.90) {
#   a <- median(x)
#   s <- sd(x)
#   n <- length(x)
#   if (n == 1){
#     s <- 0
#     n <- 2
#   }
#   error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
#   return(c(upper = a + error, med = a, lower = a - error))
# }
# 
# myCI_u <- function(x){
#   #print(length(x))
#   # x <- x[!is.na(x)]
#   bb <- med_CI(x, ci = 0.90)
#   return(as.numeric(bb[1]))
# }
# 
# myCI_l <- function(x){
#   #print(length(x))
#   # x <- x[!is.na(x)]
#   bb <- med_CI(x, ci = 0.90)
#   return(as.numeric(bb[3]))
# }
# 
# 
# ln.fxn <- function(xx, yy, pid, order.by.xx=T) {
#   #browser()
#   missing.data <- is.na(xx) | is.na(yy)
#   xx <- xx[!missing.data]
#   yy <- yy[!missing.data]
#   if(order.by.xx) {
#     ord <- order(xx)
#     xx <- xx[ord]
#     yy <- yy[ord]
#   }
#   #lines(xx, yy, col = pid[1])
# }
# 
# # =================================================================================
# 
# ## Figures 1 - Model Output
# pdf('Ch_Modelpredictions.pdf', w =10, h = 7)
# xs <- seq(0, 6, by = .1) 
# xlim <- range(xs)        
# par('ps' = 16, mfrow = c(1,1)) ## graphical parameters
# plot(0,0, type='n', xlab = 'Years since ART initiation', ylab = 'CD4 z-scores',xaxt='n', yaxt='n', bty = 'n', 
#      xlim = c(0,6000), ylim = c(-8,2), cex.main =0.9) # , main = "Suppressed viral load")
# axis(2, at = seq(-8,2, by = 2), las = 2) ## x axis
# axis(1, at = seq(0,5110, by = 730), labels = seq(0,14, by= 2) ) #,las = 2) ## x axis
# ## line for each pid
# print("Start ddply")
# #test <- ddply(modelpred, .(Ids), with, ln.fxn(years, predictions, Ids, order.by.xx=T)) 
# test <- modelpred
# test <- test[order(test$years),]                    # Order with resepct to time since HAART initiation
# dat2 <- data.matrix(test)                               # Transform dataframe into a Matrix
# length(test$lab_v) == length(dat2[,3])                 # (diff variable)
# 
# # Transform matrix into longitudinal object's class
# print("Build longitudinal object")
# dat <- as.longitudinal(dat2 , repeats = as.numeric(table(as.numeric(dat2[,3]))), unique(as.numeric(dat2[,3])))
# #is.longitudinal(dat)  
# 
# # Calculate the medians
# med <- condense.longitudinal(dat, 1, median)
# confIntu <- condense.longitudinal(dat, 1, myCI_u) #CI(dat[,33],ci = 0.95)
# confIntl <- condense.longitudinal(dat, 1, myCI_l) #CI(dat[,33],ci = 0.95)  
# tim <- get.time.repeats(dat)
# #   sp <- smooth.spline(tim$time, med, spar=0.35)
# #   lines(sp, col = col.vec[ii])
# lines(tim$time[180:length(tim$time)], rollmean(med, 180))#, col = col.vec[ii])
# # if (ii == nn){
# #   lines(tim$time[180:length(tim$time)], rollmean(confIntu,180), col =  gray(0.7), lty = 2)
# # }
# # if (ii == 1){
# #   lines(tim$time[180:length(tim$time)], rollmean(confIntl,180), col =  gray(0.7), lty = 2)
# # }
# # }
# title("CD4 z-scores medians' trajectory" ) #, outer=FALSE)
# #legend("topright", levels(testchdata$cd4a.categ), col = 1:nn, lty = 1)
# dev.off()
# 
# 
# # AIC
# anova(model1)
# anova(model1,model2)
# anova(model33,model3)   # Not comparable
# qqnorm(model1)
# qqnorm(model2)
# qqnorm(model3)
# qqnorm(model33)
# # ============================================================================
# 
