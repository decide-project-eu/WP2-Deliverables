### DECIDE Work Package 2 - Deliverable 1 - Part 1: Univariate Dynamic Linear Model
### University of Copenhagen

### Introduction to this example (according to Kristensen et al., 2010) ----

# Animal health surveillance requires a continual stream of data. It also needs a framework that allows knowledge from information previously acquired to accumulate, leading to an improved understanding of the present situation and the detection of apparent trends or unexpected changes. State-space models offer such a framework, in which relevant prior knowledge and current information are combined to detect changes in an observed process, thus allowing for a better understanding of the situation, and resulting in better-informed decisions.
# A basic model to describe a time series comprises an expected underlying value (or “true mean”), a sample error and an observation error. The fundamental assumption behind this type of model is that the true underlying mean is constant over time, but this assumption is often dubious, since values can often suffer some random fluctuations, or vary over the day or according to seasons, or just systematically increase or decrease over time.
# For that reason, a type of system is required, which allows the estimation of an observation based on the true underlying mean with its measurement or observation errors, but also of the dynamic aspects of those observations, meaning its systematic fluctuations and changes over time. Therefore, a state-space model is normally defined by two equations: an observation equation describing the data observed by a vector of parameters, and a system equation describing the dynamics of the parameter vector as a first order autocorrelation model. Components of the system equation may include, for example, fluctuations in water drinking during the day, lactation curves, disease seasonality or positive weight gain linear trends.
# For the purposes of the present project, the most commonly used type of state space model is a Dynamic Linear Model (DLM). A DLM uses a Bayesian framework to estimate the underlying parameter vector from the observed data, while taking into account any prior information available before the observations are done. Values are forecast at each time-step, based on the theoretical true mean and prior knowledge on error and variance around the system and the data, and are then updated on the next time step, being “corrected” according to each new observation, as well as to the error and variance components already mentioned.
# (Kristensen AR, Jørgensen E, Toft N. 2010. Herd Management Science II. Advanced Topics. Chapter 8.)

## In this document, we wish to walk the reader through an approach to monitor the mortality of salmon produced through aquaculture in Scotland using DLMs. Two broad methods of DLM will be applied. These two types of DLM are:
#  The univariate DLM, where a DLM is defined and applied to each variable separately
#  The multivariate DLM, where a single DLM is defined to describe all the relevant variables simultaneously

# This first example will go through the univariate DLM.

### Start by using the *load* function to load the salmon anonymized dataset:

getwd()
load("df_salmon.RData")
df <- df_salmon

# The type of data this script will be applied is to mortality data from salmon produced in Scotland (during the seawater phase) from several farms/sites. Besides the mortality variable, there are also environmental variables, e.g. sea temperature, salinity, pH...
# Each farm (in aquaculture called site) operates a single production cycle at a time (production cycle can be understood as the same as batch). Between cycles, there are fallow periods during which no fish are present on the site.

### Now, use the *source* function to source the following script that contains the DLM functions:

source('DLM functions - DECIDE deliverable_Resubmission.R')

### And we will need the following libraries:

library(dplyr)

### Make a vector called "relevant.names", which includes the following column names to be used:

relevant.names <- c("log0.mortality.rel.20")


### Create Learning and Test sets and Standardize the data ----

# The dataset is split between a Learning set (to train the model) and a Test set, to try and detect changes from the "normal" baseline obtained with the Learning set.
## Different ways can be used to create the Learning and Test sets, e.g. healthy farms go to the Learning.set and farms with reported outbreaks to the Test.set, simple division on the dataset into Learning and Test sets, being the 70:30, 75:25, 80:20 proportions the most commonly used, etc.
## We will use the first ~3/4 of time (dates) for the Learning.set (that is the 75:25 proportion). The production cycles that are cut using the 75:25 cutoff date are moved to the Learning.set if most of the production cycle in on the Learning.set already, or to the Test.set if most of the information about the production cycle is the in the Test.set.
N <- round(3*(length(unique(df[order(df[,"date"]), "date"]))/4))
sets <- get.learning.test.sets(df, relevant.names=relevant.names, hierarchical=FALSE, N=N)
Learning.set <- sets[["Learning.set"]]
Test.set <- sets[["Test.set"]]

# When variables differ too much in range and magnitude, larger numerical changes can drive the model, so it is customary to standardize them. It is not strictly necessary to standardize the data used in a **univariate**, but there is no harm in it either. So we will standardize now. 

### Make the for-loop
SDs <- c()
Means <- c()
for(name in relevant.names){
  SDs <- c(SDs, sd(na.omit(Learning.set[,name])))
  Means <- c(Means, mean(na.omit(Learning.set[,name])))
}

for(name in relevant.names){
  Mean.name <- Means[which(relevant.names == name)]
  SD.name <- SDs[which(relevant.names == name)]
  Learning.set[,name] <- (Learning.set[,name] - Mean.name)/SD.name
  Test.set[,name] <- (Test.set[,name] - Mean.name)/SD.name
}

Standarized.factors <- cbind(as.data.frame(Means), as.data.frame(SDs))
rownames(Standarized.factors) <- relevant.names


### Univariate DLM on salmon mortality ###

### Get the relevant information: mu0, C0, Gt, Ft, V and W ----

### Use the *get.mu0* (from the sourced functions file) function to make the initial parameter vector, mu0, and  *get.C0* to make the initial variance matrix, C0, and call the outputs "outmu" and "outc". Use the following settings when applying for the function:

# Data = Learning.set
# stratify.by = NA. This input will be used when you want to monitor a variable separately in different situations. As an example, we could stratify the data by production cycle/batch, by farm, or whatever other grouping that, from a biological or management point of view, could generate different data behavior groups in the monitored variable.
# time.var = 'months.since.start'. This signifies the months since salmon are stocked into the sea until they are harvested.
# expected.start.time = 0 . This determines from which value of the time variables you want to start monitoring. You could, for example, believe that in the first month of salmon production there is a lot of uncertain events to affect mortality, and therefore it is best to start monitoring from months.since.start=1 or months.since.start=2, instead of 0 (in our dataset 0 is the first month). But for this exercise, we will start from 0, since we want to start by the first information we have (which in our case if from month 0)
# relevant.names <- c('log0.mortality.rel.20')
# simple.linear = FALSE . This input is used to "tell" the function if your data should have a linear trend or not. By looking at the function file, you can see that, if the trend is linear, mu0 will be a matrix with the average value for day 1 on the first row, and the expected change for the next on the second row. Further down on the model, it will be multiplied with Gt and the updated mu with Ft, which will have the necessary structures to handle that. If the expected changes are not linear, the function will build mu0 as a matrix with mu on the first row and 1 on the second row, and Ft and Gt will also be adjusted accordingly.

# From *outmu*, extract *mu0* and from *outc" extract *C0* and call them "mu0" and "C0" respectively

outmu <- get.mu0(Data = Learning.set,
                 stratify.by=NA, 
                 time.var = 'months.since.start', 
                 expected.start.time = 0, 
                 relevant.names = c('log0.mortality.rel.20'), 
                 simple.linear = FALSE)

outc <- get.C0(Data = Learning.set,
               stratify.by=NA, 
               time.var = 'months.since.start', 
               expected.start.time = 0, 
               relevant.names = c('log0.mortality.rel.20'))

mu0 <- outmu$mu0
C0 <- outc$C0

### Use the get.Gt and get.Ft functions to create the Gt and Ft matrices, respectively. When running the two functions, set *relevant.names* = c('log0.mortality.rel.20')

### In order to run the Gt matrix, you need to first define a spline for the data series, if your DLM contains a non-linear trend component, as explained earlier. You can do that by using function get.spline.

Spline.list <- get.spline(Data = Learning.set,
                          stratify.by=NA, 
                          time.var = 'months.since.start', 
                          relevant.names = c('log0.mortality.rel.20'),
                          plot.it = FALSE)

Gt <- get.Gt(Data.A = Learning.set,
             i.A=NA,
             time.var.A='months.since.start', 
             stratify.by.A=NA, 
             Spline.list.A=Spline.list, 
             relevant.names.A=c('log0.mortality.rel.20'))

Ft <- get.Ft(relevant.names.A = c('log0.mortality.rel.20'))

### Use function *get.V* to find the observational variance, V, for *log0.mortality.rel.20*. Use the following settings:

# Data = Learning.set
# identifyer = 'nseq' (this corresponds to the production cycle ID)
# relevant.names = 'log0.mortality.rel.20'

## get.V requires de calculation of a moving average, and runs function moving.function() inside it. This function should have been loaded together with the others when sourcing the functions file.

V <- get.V(Data = Learning.set,
           identifyer='nseq',
           stratify.by=NA,
           time.var='months.since.start',
           relevant.names=c('log0.mortality.rel.20'))

### After obtaining a data-based V from the function above, it is possible to directly estimate W and update/optimize the V initially obtained with get.V, by using an Estimation Maximization (EM) algorithm. The EM function uses several functions that have been sourced from the functions file, such as getVSumElement, runDLM (we will come back to this one later), runSmoother, and the EM algorithm for a specific number of steps (runEM), since we are going to run a version of the same algorithm with an early-stopping feature that stops running when the values of V and W yielding the best model performance are found.

varcom <- runEM_earlyStopping(Data=Learning.set,
                              stratify.by=NA,
                              Spline.list=Spline.list,
                              identifyer='nseq',
                              V0=V$V_1, 
                              W0=NA,
                              C0.list=outc,
                              mu0.list=outmu,
                              no.better.limit=1,
                              time.var='months.since.start',
                              relevant.names=c('log0.mortality.rel.20'),
                              round.by=4)

## The output "varcom" should be a list containing objects V.list, W.list, mu0.list and C0.list. Those, in turn, contain V_1 (which is the EM-optimized version of the V_1 initially obtained with get.V), W_1 (the value of W estimated by the EM algorithm), mu0_1 (which is the EM-optimized version of the mu0 initially obtained with get.mu0) and C0_1 (which is the EM-optimized version of the C0 initially obtained using get.C0 and stored in outc). Here, only the optimized V and W will be used as input for your DLM on the next step, keeping the original mu0 and C0 obtained using get.mu0 and get.C0. However, all variance components given by the EM algorithm can be used when applying the DLM to the Test.set. 

## In situations for which there is not enough data, not enough time or not enough computing power to use the EM algorithm to optimize V and obtain W, it is possible to use a discount factor (delta) as the percentage of the total variance C0 that corresponds to W. As in the EM algorithm, several values of delta can be tested and optimized.
### Use the *optimize.delta* function from the functions script to find the optimal delta value for modelling mortality with the DLM. Apply the function with the following settings:
 
delta <- optimize.delta(deltas=seq(from=0.5, to=1, by=0.01), 
                        Learningset=Learning.set,
                        identifyer='nseq',
                        mu0.list=outmu,
                        C0.list=outc, 
                        V.list=V$V_1,
                        relevant.names= c('log0.mortality.rel.20'),
                        et.name="log0.mortality.rel.20",
                        Spline.list=Spline.list,
                        time.var='months.since.start',
                        stratify.by=NA)


### Validate the model ----

## Apply the optimized DLM to the mortality of each production cycle (nseq) in the Learning.set, in order to validate the model by checking whether the standardized forecast errors follow a standard normal distribution.

# Make an empty vector
extracted.all <- data.frame()

# Define the identifyer, that is, the column name that refers to the ID of each production cycle
identifyer <- 'nseq'

# Use a *for*-loop to iterate over each of the Learning.set production cycles IDs:
for(ID in unique(Learning.set[,identifyer])){
  
  # Extract the subset for the current ID
  ID.set <- subset(Learning.set, Learning.set[,identifyer] == ID)
  
  # Apply the DLM and call the result "res".
  res <- runDLM(ID.set, 
                mu0=outmu$mu0, 
                C0=outc$C0, 
                V=varcom$V.list$V_1, 
                W=varcom$W.list$W_1, 
                adjust.W=FALSE, 
                delta=NA, 
                relevant.names=c('log0.mortality.rel.20'), 
                Spline.list=Spline.list, 
                time.var='months.since.start', 
                stratify.by=NA)
  
  # Use the *extract.res* function to extract the output of the DLM (which is returned as lists) as a data frame, and call the resulting data frame "extracted".
  extracted <- extract.res(res = res, smot = NULL, relevant.names=c('log0.mortality.rel.20'))
  extracted.all <- rbind(extracted.all, extracted)
}

### Validate the DLM by assessing whether the standardized forecast errors follow a standard normal distribution, plotting it and getting the percentage outside of the 95 % CI. Use function *assess.ut* for that.
assess.ut(extracted.all)


### Apply the univariate DLM to the Test.set ----

## Now we wish to apply the DLM to the *Test.set* with the optimized parameters: mu0, C0, V, and W (V and W from the EM algorithm)

## Here we will also apply the Smoother to the DLM results (res). The retrospective smoothing analyses the data backwards, providing the best possible estimate of the true underlying mortality level given all available information prior to and after a given time step. Therefore, it gives the best estimates for mortality.

# Make an empty vector
extracted.all.test <- data.frame()

# Use a *for*-loop to iterate over each of the Test.set production cycles IDs:
for(ID in unique(Test.set[,identifyer])){
  
  # Extract the subset for the current ID
  ID.set <- subset(Test.set, Test.set[,identifyer] == ID)
  
  # Apply the DLM and call the result "res".
  res <- runDLM(ID.set, 
                mu0=outmu$mu0, 
                C0=outc$C0, 
                V=varcom$V.list$V_1, 
                W=varcom$W.list$W_1, 
                adjust.W=FALSE, 
                delta=NA, 
                relevant.names=c('log0.mortality.rel.20'), 
                Spline.list=Spline.list, 
                time.var='months.since.start', 
                stratify.by=NA)
  
  # Run the Smoother on the DLM outcomes
  smot <- runSmoother(res)
  
  # Use the *extract.res* function to extract the output of the DLM (which is returned as lists) as a data frame, and call the resulting data frame "extracted".
  extracted <- extract.res(res = res, smot = smot, relevant.names=c('log0.mortality.rel.20'))
  extracted.all.test <- rbind(extracted.all.test, extracted)
}


### Plot DLM results for all production cycles ----

par(mfrow = c(3, 1), mar = c(4, 4, 4, 2) + 0.1)

# Use a *for*-loop to iterate over each of the Test.set production cycles IDs:
## In this case, there are a lot of ID's so it may take some time. We suggest to decide on the ID's we want to see and run one by one. Run: ID=unique(Test.set[,identifyer])[353], ID=unique(Test.set[,identifyer])[71]
for(ID in unique(Test.set[,identifyer])){
  
  ID.set <- subset(Test.set, Test.set[,identifyer] == ID)
  i <- which(Test.set[,identifyer] == ID)
  mt <- extracted.all.test[i, paste("mt_", relevant.names, sep="")] 
  ft <- extracted.all.test[i, paste("ft_", relevant.names, sep="")]
  Ct <- extracted.all.test[i, paste("Ct_", relevant.names, sep="")]      
  Qt <- extracted.all.test[i, paste("Qt_", relevant.names, sep="")]      
  mts <- extracted.all.test[i, paste("mts_", relevant.names, sep="")]     
  Cts <- extracted.all.test[i, paste("Cts_", relevant.names, sep="")]     
  
  min <- min(ft-1.96*sqrt(Qt), na.rm=T)
  max <- max(ft+1.96*sqrt(Qt), na.rm=T)
  
  #mt - filtered mean
  plot(ID.set[,relevant.names], type="l", xlab="month", ylab="log.mortality", 
       main=c(paste("Filtered mean - nseq:", ID, sep=" ")), 
       ylim=c(min,max), lwd=2)                           # Observations
 
  low.limit.mt <- mt-1.96*sqrt(Ct)
  high.limit.mt <- mt+1.96*sqrt(Ct)
  lines(mt, col='darkgreen', lwd=2)                      # Filtered mean
  lines(low.limit.mt, col='darkgreen', lty = "dashed")   # Filtered variance
  lines(high.limit.mt, col='darkgreen', lty = "dashed")  # Filtered variance
  
  legend(x = "bottomright",                              # Position
         legend = c("Obs", "mt", "95% CI"),              # Legend texts
         col = c("black", "darkgreen", "darkgreen"),     # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=1)
  
  #ft - forecasts (if above the credible intervals it can be considered that mortality is significantly higher than expected -> warning)
  plot(ID.set[,relevant.names], type="l", xlab="month", ylab="log.mortality", 
       main=paste("Forecasts - nseq:", ID, sep=" "), 
       ylim=c(min,max), lwd=2)                           # Observations

  low.limit.ft <- ft-1.96*sqrt(Qt)
  high.limit.ft <- ft+1.96*sqrt(Qt)
  lines(ft, col='red', lwd=2)                            # Forecasts
  lines(low.limit.ft, col='red', lty = "dashed")         # Forecast variance
  lines(high.limit.ft, col='red', lty = "dashed")        # Forecast variance
  
  legend(x = "bottomright",                              # Position
         legend = c("Obs", "ft", "95% CI"),              # Legend texts
         col = c("black", "red", "red"),                 # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=1)
  
  #mts - smoothed mean
  plot(ID.set[,relevant.names], type="l", xlab="month", ylab="log.mortality", 
       main=paste("Smoothed mean - nseq:", ID, sep=" "), 
       ylim=c(min,max), lwd=2)                           # Observations
  
  low.limit.mts <- mts-1.96*sqrt(Cts)
  high.limit.mts <- mts+1.96*sqrt(Cts)
  lines(mts, col='blue', lwd=2)                          # Smooth mean
  lines(low.limit.mts, col='blue', lty = "dashed")       # Smooth variance
  lines(high.limit.mts, col='blue', lty = "dashed")      # Smooth variance
  
  legend(x = "bottomright",                              # Position
         legend = c("Obs", "mts", "95% CI"),             # Legend texts
         col = c("black", "blue", "blue"),               # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=1) 
  
  ### See real values ###
  
  # - remove standardization from mortality
  Mean.name <- Standarized.factors[relevant.names, "Means"]
  SD.name <- Standarized.factors[relevant.names, "SDs"]
  
  Obs_no_s <- subset(Test.set[,relevant.names], Test.set[,identifyer] == ID)
  mt_no_s <- (mt * SD.name) + Mean.name
  ft_no_s <- (ft * SD.name) + Mean.name
  mts_no_s <- (mts * SD.name) + Mean.name
  
  ## - For the CI's we have to remove the standardization on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt_no_s <- (low.limit.mt * SD.name) + Mean.name
  high.limit.mt_no_s <- (high.limit.mt * SD.name) + Mean.name
  low.limit.ft_no_s <- (low.limit.ft * SD.name) + Mean.name
  high.limit.ft_no_s <- (high.limit.ft * SD.name) + Mean.name
  low.limit.mts_no_s <- (low.limit.mts * SD.name) + Mean.name
  high.limit.mts_no_s <- (high.limit.mts * SD.name) + Mean.name
  
  # - remove log transformation from mortality
  Obs <- subset(Test.set[,"mortality.rel.20"], Test.set[,identifyer] == ID)
  mt <- exp(1)^mt_no_s - 1/20000
  ft <- exp(1)^ft_no_s - 1/20000
  mts <- exp(1)^mts_no_s - 1/20000
  
  ## - For the CI's we have again to remove log transformation on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
  high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
  low.limit.ft <- exp(1)^low.limit.ft_no_s - 1/20000
  high.limit.ft <- exp(1)^high.limit.ft_no_s - 1/20000
  low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
  high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000
  
  # - plot with the real values
  min <- min(c(low.limit.ft, Obs), na.rm=T)
  max <- max(c(high.limit.ft, Obs), na.rm=T)
  
  #mt - filtered mean
  plot(Obs, type="l", xlab="month", ylab="mortality",
       main=paste("Filtered mean - nseq:", ID, sep=" "),
       ylim=c(min,max), lwd=2)                           # Observations
  
  lines(mt, col='darkgreen', lwd=2)                      # Filtered mean
  lines(low.limit.mt, col='darkgreen', lty = "dashed")   # Filtered variance
  lines(high.limit.mt, col='darkgreen', lty = "dashed")  # Filtered variance
  
  legend(x = "topright",                                 # Position
         legend = c("Obs", "mt", "95% CI"),              # Legend texts
         col = c("black", "darkgreen", "darkgreen"),     # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=1)
  
  #ft - forecasts (if above the credible intervals it can be considered that mortality is significantly higher than expected -> warning)
  plot(Obs, type="l", xlab="month", ylab="mortality",
       main=paste("Forecasts - nseq:", ID, sep=" "), 
       ylim=c(min,max), lwd=2)                           # Observations
  
  lines(ft, col='red', lwd=2)                            # Forecasts
  lines(low.limit.ft, col='red', lty = "dashed")         # Forecast variance
  lines(high.limit.ft, col='red', lty = "dashed")        # Forecast variance
  
  legend(x = "topright",                                 # Position
         legend = c("Obs", "ft", "95% CI"),              # Legend texts
         col = c("black", "red", "red"),                 # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=1)
  
  #mts - smoothed mean
  plot(Obs, type="l", xlab="month", ylab="mortality",
       main=paste("Smoothed mean - nseq:", ID, sep=" "),
       ylim=c(min,max), lwd=2)                           # Observations
  
  lines(mts, col='blue', lwd=2)                          # Smooth mean
  lines(low.limit.mts, col='blue', lty = "dashed")       # Smooth variance
  lines(high.limit.mts, col='blue', lty = "dashed")      # Smooth variance
  
  legend(x = "topright",                                 # Position
         legend = c("Obs", "mts", "95% CI"),             # Legend texts
         col = c("black", "blue", "blue"),               # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=1)
}
