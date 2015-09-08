rm(list = ls())

# libraries
library(nlme)

######################## FUNCTIONS   ########################################

# Select a subgroup of people with TOFO >= 8 years
useSubset <- T
subselect <- function(addata,n){ 
  num.to.do <- n ## to save computing time when playing with analyses
  pids.to.run1 <- sample(unique(addata$patient), num.to.do) 
  addata <- addata[addata$patient %in% pids.to.run1,] 
  return(addata)
}

# Select a subgroup of people
useSubset <- T
ltofuselect <- function(addata,n){ 
  num.to.do <- n ## to save computing time when playing with analyses
  pids.to.run1 <- sample(unique(addata$patient[addata$diff > 2920]), num.to.do) 
  addata <- addata[addata$patient %in% pids.to.run1,] 
  return(addata)
}

# Likelihood for the Asymptotic model
negloglkhd <- function(parms,y) {   ## Rosenbrock Banana function
  yy <- y$lab_v
  tim <- y$diff
  # Define the parameters
  K = parms[1]
  Q = parms[2]
  r = parms[3]
  s = parms[4]
  z0 = parms[5]
  
  if( K == z0){
    print( "negloglkhd cannot be evaluated with given parameters: K = z0" )
  }  
  else if (Q == 1){
    print( "negloglkhd cannot be evaluated with given parameter: Q = 1" )
  }  
  else {
    fx <- yy - ( (K*z0*(Q-1)/Q*(K-z0))*(exp((r - s)*tim) + exp(s*tim)/(Q - 1))/(1 + exp(r*tim)*z0 /(K - z0)) )
    #fx <- yy - ( (K/Q) * ( (1 + exp(-s*tim)/(Q - 1) )/(1 + exp(-r*tim)*z0 /(K - z0)) ) ) 
    return( (sum(fx))^2 )
  }  
}

indfit <- function(ydat,n){
  nofvis <- length(ydat$lab_v)          # num of obs
  mat <- t(sapply(1:n, indoptim, ydat ))
  out <- mat[match(min(mat[,1]), mat),] # extract the output with the minimun RSS
  #print(out)
  # Procede with the derivative parameter estimation.
  tryCatch(res <- optim(out[2:6], negloglkhd, y = ydat,  method = "BFGS", control = c(trace = FALSE) ), 
           error = function(e) res = c("The optimization Failed") )
  
  # Results if derivative optm fails
  if ( class(res) == "character" ){ 
    sol <- out[1] / nofvis              # scale RSS
    resout <- c(sol, out[2:6],1,nofvis)
  }
  # Results from derivative optimization
  else if ( class(res) ==  "list" ){ 
    sol <- res$value / nofvis
    resout <- c(sol, res$par,2,nofvis)
  }
  return ( resout )
}

getmedpar <- function(dat,n,m){
  matn <- data.frame(RSS= double(),K =double(),Q= double(),r = double(), s = double(),z0 = double(),Optm = double(), novis = double())
  for (i in 1:m) { matn[i,] <- indfit(dat,n) } # { mat3[i,] <- indfit(dat,30)}
  best <- matn[median(matn$RSS),]
  return (best)
}


indoptim <- function(i,optdata){
  # Does one single SANN optimization
  out <- optim(c(5,3.3,0.2,0.1,2), negloglkhd, y = optdata,  method = "SANN", control = c(trace = FALSE) )
  #lower = c(-10,-10,-5), upper = c(5,5,5) )
  return (c(out$value, out$par) )  
}

divind <- function(data,n,m){
  # data need to be a grouped dataset otherwise the function fails.
  level <- length(getGroupsFormula(data, asList = TRUE))
  groups <- getGroups(data, level = level)[drop = TRUE]
  res <- lapply(split(data, groups), getmedpar, n,m) 
  test <- as.data.frame(do.call(rbind, res))         # Define a dataframe with the output from the Optimization
  test <- cbind(patient = rownames(test), test)      # Add patient IDs to the dataframe
  names(test) <- c("patient","RSS","K","Q","r", "s","z0", "Optm", "novis")     # Change columns' names
  return(test)
}

########################################################################

# Read Data
chdata <- read.csv("/home/evaliliane/Documents/PhD/Codes/NewData/CD4Cat_Children2015-04-20.csv")

# Removing row with TOFU > 15
chdata <- chdata[chdata$diff < 5476,]
dat <- chdata[,c(2,44,35,49)]
dat <- dat[!(is.na(dat$cd4a.categ)),]
dat <- dat[!(is.na(dat$lab_v)),]
dat <- dat[!(is.na(dat$patient)),]
dat <- dat[!(is.na(dat$diff)),]

# Order the dataframe and run the individual fits with modified function.
mydata <- groupedData(lab_v ~ diff | patient, data = dat ,order.groups=F)  

##########################t## MLE ####################################################

# Test small sample
optdata <- subselect(mydata,10)

setwd("/home/evaliliane/Documents/PhD/Codes")

# Run full children dataset and save output
system.time(indout <- divind(optdata,5,5))

system.time(indout2 <- divind(optdata,5,5))
#indout <- divind(mydata,10)


write.csv(indout, file = "IndividualData10-New.csv")
