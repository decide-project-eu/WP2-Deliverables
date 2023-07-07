### DECIDE Work Package 2 - Deliverable 1 - Part 2: Multivariate Dynamic Linear Model
### University of Copenhagen
### 2023-06-26

### Introduction to this example (according to Kristensen et al., 2010) ###

# Animal health surveillance requires a continual stream of data. It also needs a framework that allows knowledge from information previously acquired to accumulate, leading to an improved understanding of the present situation and the detection of apparent trends or unexpected changes. State-space models offer such a framework, in which relevant prior knowledge and current information are combined to detect changes in an observed process, thus allowing for a better understanding of the situation, and resulting in better-informed decisions.
# A basic model to describe a time series comprises an expected underlying value (or “true mean”), a sample error and an observation error. The fundamental assumption behind this type of model is that the true underlying mean is constant over time, but this assumption is often dubious, since values can often suffer some random fluctuations, or vary over the day or according to seasons, or just systematically increase or decrease over time.
# For that reason, a type of system is required, which allows the estimation of an observation based on the true underlying mean with its measurement or observation errors, but also of the dynamic aspects of those observations, meaning its systematic fluctuations and changes over time. Therefore, a state-space model is normally defined by two equations: an observation equation describing the data observed by a vector of parameters, and a system equation describing the dynamics of the parameter vector as a first order autocorrelation model. Components of the system equation may include, for example, fluctuations in water drinking during the day, lactation curves, disease seasonality or positive weight gain linear trends.
# For the purposes of the present project, the most commonly used type of state space model is a Dynamic Linear Model (DLM). A DLM uses a Bayesian framework to estimate the underlying parameter vector from the observed data, while taking into account any prior information available before the observations are done. Values are forecast at each time-step, based on the theoretical true mean and prior knowledge on error and variance around the system and the data, and are then updated on the next time step, being “corrected” according to each new observation, as well as to the error and variance components already mentioned.
# (Kristensen AR, Jørgensen E, Toft N. 2010. Herd Management Science II. Advanced Topics. Chapter 8.)

## In this document, we wish to walk the reader through an approach to detect mastitis from the standardized forecast errors from DLMs. Two broad methods of DLM will be applied. These two types of DLM are:
#  The univariate DLM, where a DLM is defined and applied to each variable separately
#  The multivariate DLM, where a single DLM is defined to describe all the relevant variables simultaneously

# The following code will apply a multivariate DLM to the individual cow-parities.

### Start by using the *source* function to source the following script:

getwd()
source('DLM functions - DECIDE deliverable.R')

# The type of data this script would be applied to is milk yield from dairy cows in a farm where mastitis is detected every now and then. Besides the variables containing yield per se, there are also variables describing milk composition, SCC and blood presence.
# The observational units are cow-parities, meaning that she same cow will appear as two different individuals  at two different lactations. This allows for the separation of the data per lactation group, if necessary, and also isolates the developing production curve a cow has in an individual lactation. Milk yield measurements are presented on a daily basis, so the time unit is days, and the time-step unit is days in milking (dim). Cows on the first lactation are defined as "primiparous" cows, against "multiparous" cows from other lacatations. This information is contained in a variable named "ParityGroup".

# The dataset is split between healthy cows (to train the model) and sick cows, to try and detect changes from the "normal" baseline obtained with the healthy cow data. In this code, those files are called "healthy.set" and "sick.set", respectively.

# For the purposes of this report, we will only consider "Primiparous" cows. 
# Use the *subset* function the extract the subset of *healthy.set* where *ParityGroup* == "Primiparous". Call the resulting data frame "healthy.set" (i.e. replace the existing data frame with the subset you extract). do the same for the "sick.set".

healthy.set <- subset(healthy.set, healthy.set$ParityGroup == "Primiparous")
sick.set <- subset(sick.set, sick.set$ParityGroup == "Primiparous")

# When variables differ too much in range and magnitude, larger numerical changes can drive the model, so it is customary to standardize them. It is not strictly necessary to standardize the data used in a **univariate**, but there is no harm in it either. So we will standardize all variables now. 

# Make a vector called "names", which includes the following column names:

# "yield"
# "conductivity" 
# "fat"
# "protein"
# "lactose" 
# "SCC"
# "blood"

names <- c("yield",
           "conductivity",
           "fat",
           "protein",
           "lactose", 
           "scc",
           "blood")

# Make the for-loop
for(name in names){
  # Get mean and sd from the healthy.set
  Mean <- mean(healthy.set[,name])
  SD <- sd(healthy.set[,name])
  
  #Standardize the data in the healthy.set
  healthy.set[,name] <- (healthy.set[,name] - Mean)/SD
  
  #Standardize the data in the sick.set
  sick.set[,name] <- (sick.set[,name] - Mean)/SD
}

# Use the *subset* function to extract the subset of ***healthy.set*** where *dim* <= 5, and call resulting data frame "healthy_start.set".

healthy_start.set <- subset(healthy.set, healthy.set$dim <= 5)

# Use the *unique* function the get the unique names of the cow-parities in the *healthy.set* and the *sick.set*, and call the resulting vectors "IDs.healthy" and "IDs.sick", respectively. You will use these vectors later. 

IDs.healthy <- unique(healthy.set$Cow.Parity)
IDs.sick <- unique(sick.set$Cow.Parity)

# Create a new data frame called "meta.data.sick", consisting of the following columns from the "sick.set":

# "Cow.Parity"   
# "cow"          
# "lactno"       
# "ParityGroup"  
# "dim" 
# "mastitis"

meta.data.sick <- sick.set[, c("Cow.Parity", "cow", "lactno", "ParityGroup", "dim", "mastitis")]

#############################################################################

## Multivariate DLM on milk yield ##

# Use the *get.mu0_C0* (from the sourced functions file) function to make the initial parameter vector, mu0, and initial variance matrix, C0, for milk yield, and call the output "out". Use the following settings when applying for the function:

# Data = healthy_start.set
# stratify.by = NA . This input will be used when you want to monitor a variable separately in different situations. As an example, we could stratify the data by lactation phase, or by herd, or whatever other grouping that, from a biological or management point of view, could generate different data behaviour groups in the monitored variable. 
# time.var = 'dim'
# expected.start.time = 1 . This determines from which value of the time variables you want to start monitoring. You could, for example, believe that in the first week of lactation there is a lot of uncertain events to affect lactation, and therefore it is best to start monitoring from dim=5 or dim=7, instead of 1. But for this exercise, we will start from 1, since we are only looking at the frist 5 days anyway.
# relevant.names <- c('yield')
# simple.linear = FALSE . This input is used to "tell" the function if your data should have a linear trend or not. By looking at the function file, you can see that, if the trend is linear, mu0 will be a matrix with the average value for day 1 on the first row, and the expected change for the next on the second row. Further down on the model, it will be multiplied with Gt and the updated mu with Ft, which will have the necessary structures to handle that. If the expoected changes are not linear, the function will build mu0 as a matrix with mu on the first row and 1 on the second row, and Ft and Gt will also be adjusted accordingly.

# From *outmu*, extract *mu0* and from *outc" extract *C0* and call them "mu0" and "C0" respectively

outmu <- get.mu0(Data = healthy_start.set,
                        stratify.by=NA, 
                        time.var = 'dim', 
                        expected.start.time = 1, 
                        relevant.names = names, 
                        simple.linear = FALSE)

outc <- get.C0(Data = healthy_start.set,
               stratify.by=NA, 
               time.var = 'dim', 
               expected.start.time = 1, 
               relevant.names = names)

mu0 <- outmu$mu0
C0 <- outc$C0

# Use the get.Gt and get.Ft functions to create the Gt and Ft matrices, respectively. When running the two functions, set *relevant.names* = c('yield')

# In order to run the Gt matrix, you need to first define a spline for the data series, if your DLM contains a non-linear trend component, as explained earlier. You can do that by using function get.spline.

Spline.list <- get.spline(Data = healthy.set,
                          stratify.by=NA, 
                          time.var = 'dim', 
                          relevant.names = names,
                          plot.it = FALSE)


Gt <- get.Gt(Data.A=healthy.set,
             i.A=NA,
             time.var.A='dim', 
             stratify.by.A=NA, 
             Spline.list.A=Spline.list, 
             relevant.names.A=names)

Ft <- get.Ft(relevant.names.A = names)

# Use function *get.V* to find the observational variance, V, for *yield*. Use the following settings:

# Data = healthy.set
# ID.var = 'Cow.Parity'
# relevant.names = 'yield'

# get.V requires de calculation of a moving average, and runs function moving.function() inside it. This function should have been loaded together with the others when sourcing the functions file.

V <- get.V(Data = healthy.set, 
           identifyer='Cow.Parity', 
           stratify.by=NA, 
           time.var='dim', 
           relevant.names=names)

# After obtaining a data-based V from the function above, it is possible to directly estimate W and update/optimize the V initially obtained with get V, by using an Estimation Maximization (EM) algorithm. The EM function uses several functions that have been sourced from the functions file, such as getVSumElement, runDLM (we will come back to this one later), runSmoother, and the EM algorythm for a specific number of steps (runEM), since we are going to run a version of the same algorythm with an early-stopping feature that stops running when the values of V and W yielding the best model performace are found.

varcom <- runEM_earlyStopping(Data=healthy.set, 
                              stratify.by=NA, 
                              Spline.list=Spline.list, 
                              identifyer='Cow.Parity', 
                              V0=V$V_1, 
                              W0=NA, 
                              C0.list=outc, 
                              mu0.list=outmu, 
                              no.better.limit=1, 
                              time.var='dim', 
                              relevant.names=names, 
                              round.by=4)

# The output "varcom" should be a list containing objects V.list, W.list, muo.list and C0.list. Those, in turn, contain V_1 (which is the EM-optimized version of the V_1 initially obtained with get.V), W_1 (the value of W estimated by the EM algorithm), mu0_1 (the same mu0 obtained earlier by get.mu0 and stored in mu0) and C0_1 (the same C0 obtained earlier using get.C0 and stored in outc). Those should be the variance components used as input for your DLM on the next step.

# In situations for which there is not enough data, not enough time or not enough computing power to use the EM algorithm to optimize V and obtain W, it is possible to use a discount factor (delta) as the percentage of the total variance C that corresponds to W. As in the EM algorithm, several values of delta can be tested and optimized.
  
# Use the *optimize.delta* function from the functions script to find the optimal delta value for modelling yield with the DLM. Apply the function with the following settings:

delta <- optimize.delta(deltas=seq(from=0.8, to=1, by=0.01), 
                        healthy.set, 
                        identifyer= 'Cow.Parity', 
                        mu0.list=outmu, 
                        C0.list=outc,  
                        V.list=V$V_1, 
                        relevant.names= names, 
                        Spline.list=Spline.list, 
                        time.var='dim', 
                        stratify.by=NA)

# Apply the optimized DLM to the yield of each cow-parity in the healthy set, in order to estimate the parameters we need for the control chart. 

# Make an empty vector for the standardized forecast errors of all the healthy cow-parities
ut.all <- c() 
extracted.all <- data.frame()

# Use a *for*-loop to iterate over each of the healthy IDs; for ID in IDs.healthy, do the following:

for(ID in IDs.healthy[c(1:3)]){
  
  # Extract the subset for the current ID
  ID.set <- subset(healthy.set, healthy.set$Cow.Parity == ID)
  
  # Apply the DLM and call the result "res"
  
  res <- runDLM(ID.set, 
                mu0=varcom$mu0.list$mu0_1, 
                C0=varcom$C0.list$C0_1, 
                V=varcom$V.list$V_1, 
                W=varcom$W.list$W_1, 
                adjust.W=FALSE, 
                delta=NA, 
                relevant.names=names, 
                Spline.list=Spline.list, 
                time.var='dim', stratify.by=NA)
 
  # Use the *extract.res* function to extract the output of the DLM (which is returned as lists) as a data frame, and call the resulting data frame "extracted".
  
  extracted <- extract.res(res = res, relevant.names=names)
  extracted.all <- rbind(extracted.all, extracted)
  
  # From *extracted*, get the standardized forecast errors from the column called "ut_yield" and call it "ut".
  ut <- data.frame(extracted[,grep(pattern = 'ut_', x = colnames(extracted))])
  colnames(ut) <- colnames(extracted)[grep(pattern = 'ut_', x = colnames(extracted))]
  
  # Add the standardized forecast errors to the ut.all vector
  ut.all <- rbind(ut.all, ut)
 
 
}

# Validate the DLM by assessing whether the standardized forecast errors follow a standard normal distribution, plotting it and getting the percentage outside of the 95 % CI. Use function *assess.ut* for that.

assess.ut(extracted.all)

# Now we wish to apply the DLM with the optimized parameters for yield (mu0, C0, V, and delta) to the *sick.set*. 

ut.all.sick <- c()
extracted.all.sick <- data.frame()

for(ID in IDs.sick[c(1:4)]){
  
  ID.set <- subset(sick.set, sick.set$Cow.Parity == ID)
  
  res <- runDLM(ID.set, 
                mu0=varcom$mu0.list$mu0_1, 
                C0=varcom$C0.list$C0_1, 
                V=varcom$V.list$V_1, 
                W=varcom$W.list$W_1, 
                adjust.W=FALSE, 
                delta=NA, 
                relevant.names=names, 
                Spline.list=Spline.list, 
                time.var='dim', 
                stratify.by=NA)
  
  extracted <- extract.res(res = res, relevant.names=names)
  extracted.all.sick <- rbind(extracted.all.sick, extracted)
  
  ut <- data.frame(extracted[,grep(pattern = 'ut_', x = colnames(extracted))])
  colnames(ut) <- colnames(extracted)[grep(pattern = 'ut_', x = colnames(extracted))]
  
  ut.all <- rbind(ut.all, ut)
  
}

# We will not test the performances of each variable on their own again, because in this specific case, the performances would be comparable to what we saw before in the univariate example. Instead, we will save the data that we can extract from the output of the DLM.
