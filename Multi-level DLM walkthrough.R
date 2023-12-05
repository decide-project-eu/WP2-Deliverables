#### Univariate Multi-level DLM - DECIDE Deliverable 2.2
## Example: Scottish Salmon study (mortality)

#libraries
library(dplyr)

#get anonymized dataset
setwd("H:/PhD KU/DECIDE/Deliverable D2.2/R studio")
load("df_salmon.RData")
df <- df_salmon


### Create Learning and Test sets ----

source("~/DECIDE/Deliverable D2.2/R studio/Functions - Multi-level DLM.R")

#get Learning and Test sets
#use the first ~3/4 (that will be 70%) of time (dates) for Learning.set
N <- round(3*(length(unique(df[order(df[,"date"]), "date"]))/4))
sets <- get.learning.test.sets(df, relevant.names="log0.mortality.rel.20", hierarchical=TRUE, N=N)
Learning.set.i <- sets[["Learning.set"]] #19198 rows
Test.set.i <- sets[["Test.set"]] #7796 rows


### Re-arrange datasets to work with this Multi-level approach - one variable: Mortality log transformed as log(x+0.00005) (log0.mortality.rel.20)
#Learning.set
col.needed <- Learning.set.i[,c(1,2,3,11)] #log0.mortality.rel.20
col.needed <- col.needed %>% # - order by region
  arrange(local.authority)
col.needed[,"site_region"] <- paste(col.needed$site, col.needed$local.authority, sep="_")
col.needed <- col.needed[,c(1,4,5)]

# get re-arranged Learning.set
Learning.set <- as.data.frame(matrix(NA, length(unique(col.needed$date)), length(unique(col.needed$site_region))+1))
colnames(Learning.set) <- c("date", unique(col.needed$site_region))
Learning.set$date <- unique(col.needed$date)
Learning.set <- Learning.set %>% # - order by date
  arrange(date)
for (i in 1:dim(col.needed)[1]){ 
  Learning.set[which(Learning.set$date == col.needed[i,"date"]), col.needed[i,"site_region"]] <- col.needed[i,"log0.mortality.rel.20"]
}

#Test.set
col.needed <- Test.set.i[,c(1,2,3,11)] #log0.mortality.rel.20
col.needed <- col.needed %>% # - order by region
  arrange(local.authority)
col.needed[,"site_region"] <- paste(col.needed$site, col.needed$local.authority, sep="_")
col.needed <- col.needed[,c(1,4,5)]

# get re-arranged Test.set
Test.set <- as.data.frame(matrix(NA, length(unique(col.needed$date)), length(unique(col.needed$site_region))+1))
colnames(Test.set) <- c("date", unique(col.needed$site_region))
Test.set$date <- unique(col.needed$date)
Test.set <- Test.set %>% # - order by date
  arrange(date)
for (i in 1:dim(col.needed)[1]){ 
  Test.set[which(Test.set$date == col.needed[i,"date"]), col.needed[i,"site_region"]] <- col.needed[i,"log0.mortality.rel.20"]
}


### Create a table with information about months since start for each production cycle
#Learning.set
months.learning.set <- data.frame(matrix(NA, nrow = dim(Learning.set)[1], ncol = dim(Learning.set)[2]))
months.learning.set[,1] <- Learning.set$date
colnames(months.learning.set) <- colnames(Learning.set)

for (farm in unique(Learning.set.i$site)){ 
  df.farm <- subset(Learning.set.i, Learning.set.i$site == farm)
  
  for(date in unique(df.farm$date)){
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      months.learning.set[which(months.learning.set$date==date), 
                          paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      months.learning.set[which(months.learning.set$date==date), 
                          paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "months.since.start"]
    }
  }
}

#Test.set
months.test.set <- data.frame(matrix(NA, nrow = dim(Test.set)[1], ncol = dim(Test.set)[2]))
months.test.set[,1] <- Test.set$date
colnames(months.test.set) <- colnames(Test.set)

for (farm in unique(Test.set.i$site)){ 
  df.farm <- subset(Test.set.i, Test.set.i$site == farm)
  
  for(date in unique(df.farm$date)){ 
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      months.test.set[which(months.test.set$date==date), 
                      paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      months.test.set[which(months.test.set$date==date), 
                      paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "months.since.start"]
    }
  }
}

### Create a table with information about the number of the production cycle
#Learning.set
nseq.learning.set <- data.frame(matrix(NA, nrow = dim(Learning.set)[1], ncol = dim(Learning.set)[2]))
nseq.learning.set[,1] <- Learning.set$date
colnames(nseq.learning.set) <- colnames(Learning.set)

for (farm in unique(Learning.set.i$site)){ 
  df.farm <- subset(Learning.set.i, Learning.set.i$site == farm)
  
  for(date in unique(df.farm$date)){ 
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      nseq.learning.set[which(nseq.learning.set$date==date), 
                        paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      nseq.learning.set[which(nseq.learning.set$date==date), 
                        paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "nseq"]
    }
  }
}

#Test.set
nseq.test.set <- data.frame(matrix(NA, nrow = dim(Test.set)[1], ncol = dim(Test.set)[2]))
nseq.test.set[,1] <- Test.set$date
colnames(nseq.test.set) <- colnames(Test.set)

for (farm in unique(Test.set.i$site)){ 
  df.farm <- subset(Test.set.i, Test.set.i$site == farm)
  
  for(date in unique(df.farm$date)){ 
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      nseq.test.set[which(nseq.test.set$date==date), 
                    paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      nseq.test.set[which(nseq.test.set$date==date), 
                    paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "nseq"]
    }
  }
}



### Get Country, Region and Farms variances and discount factors, mu0 and C0  ----

# Make sure you avoid scientific notation in your output - it will make things easier!
options(scipen=999)

# Standardize the learning set data (force the data to be standard normal distributed)
all.values <- unlist(Learning.set[,c(2:ncol(Learning.set))])
Mean <- mean(na.omit(all.values))
SD <- sd(na.omit(all.values))
Standarized.factors <- c(Mean=Mean, SD=SD)

Learning.set.stand <- data.frame(matrix(NA, nrow = dim(Learning.set)[1], ncol = dim(Learning.set)[2]))
Learning.set.stand[,1] <- Learning.set$date
colnames(Learning.set.stand) <- colnames(Learning.set)

Learning.set.stand[1:dim(Learning.set.stand)[1], 2:dim(Learning.set.stand)[2]] <- (unlist(Learning.set[,c(2:ncol(Learning.set))]) - Mean)/SD

# Get the initial parameter vector (mu0), the prior variance (C0) and the Spline functions
mu <- get.m0(D=Learning.set.stand, D.months=months.learning.set, expected.start.time=0, regions=c("A", "B", "C", "D", "E"))
mu0 <- mu$mu0
Spline.list  <- mu$Spline.list

C <- get.C0(D=Learning.set.stand, D.months=months.learning.set, expected.start.time=0, regions=c("A", "B", "C", "D", "E"))
C0 <- C$C0

# Run the optimum function with DLM to get 3 V's and 3 deltas (Country, Region and Farm) - can take some time to run!!
Start.time <- Sys.time()

Vs.deltas <- estimateDiscountModel(priorVs=c(3, 2, 1), priorDeltas=c(0.99, 0.99, 0.99), 
                                  D=Learning.set.stand, D.months=months.learning.set, 
                                  mu0, C0, Spline.list, regions=c("A", "B", "C", "D", "E"))

print(Sys.time()-Start.time)

# - Save the 3 V's and 3 deltas so that you don't have to run it every time
#saveRDS(object = Vs.deltas, file = paste(model.dir, '/Vs.deltas.RDS', sep=''))

# - Get the transformed back values
Vs.deltas.t <- transformResults(Vs.deltas, vars=3)

# - Get the optimum V's and deltas
countryVar = Vs.deltas.t[1]
regionVar = Vs.deltas.t[2]
farmVar = Vs.deltas.t[3]
deltas = Vs.deltas.t[4:6]



### Apply Univariate Multi-level DLM ----

# Standardize the test set data (force the data to be standard normal distributed)
Test.set.stand <- data.frame(matrix(NA, nrow = dim(Test.set)[1], ncol = dim(Test.set)[2]))
Test.set.stand[,1] <- Test.set$date
colnames(Test.set.stand) <- colnames(Test.set)

Test.set.stand[1:dim(Test.set.stand)[1], 2:dim(Test.set.stand)[2]] <- (unlist(Test.set[,c(2:ncol(Test.set))]) - Standarized.factors["Mean"])/Standarized.factors["SD"]

# Run the DLM with the optimum V's and deltas
out <- runDiscountDLM(D=Test.set.stand, D.months=months.test.set, mu0, C0, countryVar, regionVar, 
                      farmVar, deltas, Spline.list, regions=c("A", "B", "C", "D", "E"))

# Run the Smoother on the DLM outcomes
smot <- runSmoother(out)

# Get the results for out and smot
results <- extract.res(out, smot, D=Test.set.stand, H.w.country=NULL, regions=c("A", "B", "C", "D", "E"))

# Plot the standardized forecast errors - see if follow a normal distribution 
ut.res <- results[["ut"]]
ut <- unlist(ut.res[,c(8:ncol(ut.res))])  
hist(ut)



### Plot Univariate Multi-level DLM results ----

# D is now the standardize test set
D <- Test.set.stand


## Country ----
# - Plot with log transformed and standardized values
min <- min(results$mt$Country-1.96*sqrt(results$Ct$Country), na.rm=T)
max <- max(results$mt$Country+1.96*sqrt(results$Ct$Country), na.rm=T)

#mt - darkgreen
low.limit.mt <- (results$mt$Country)-1.96*sqrt(results$Ct$Country)
high.limit.mt <- (results$mt$Country)+1.96*sqrt(results$Ct$Country)


plot(results$mt$Country~results$mt$date, type="l", xlab="date", ylab="log.mortality", 
     main="Filtered mean - Country", ylim=c(min,max), lwd=2, col='darkgreen') 
lines(low.limit.mt~results$mt$date, col='darkgreen', lty = "dashed")
lines(high.limit.mt~results$mt$date, col='darkgreen', lty = "dashed") 
legend(x = "bottomright",                    
       legend = c("mt", "95% CI"),  
       col = c("darkgreen", "darkgreen"),   
       lwd = 1,                             
       lty = c(1,2),
       cex=0.7)


#mts (smoother) - blue
low.limit.mts <- results$mts$Country-1.96*sqrt(results$Cts$Country)
high.limit.mts <- results$mts$Country+1.96*sqrt(results$Cts$Country)

plot(results$mts$Country~results$mts$date, type="l", xlab="date", ylab="log.mortality", 
     main="Smoothed mean - Country", ylim=c(min,max), lwd=2, col='blue') 
lines(low.limit.mts~results$mts$date, col='blue', lty = "dashed")
lines(high.limit.mts~results$mts$date, col='blue', lty = "dashed") 
legend(x = "bottomright",                    
       legend = c("mts", "95% CI"),  
       col = c("blue", "blue"),   
       lwd = 1,                             
       lty = c(1,2),
       cex=0.7)


# - remove standardization from mortality
mt_no_s <- (results$mt$Country * Standarized.factors["SD"]) + Standarized.factors["Mean"]
mts_no_s <- (results$mts$Country * Standarized.factors["SD"]) + Standarized.factors["Mean"]

## - For the CI's we have to remove the standardization on the all CI calculation, 
## - not individually on Ct, Qt and Cts
low.limit.mt_no_s <- (low.limit.mt * Standarized.factors["SD"]) + Standarized.factors["Mean"]
high.limit.mt_no_s <- (high.limit.mt * Standarized.factors["SD"]) + Standarized.factors["Mean"]
low.limit.mts_no_s <- (low.limit.mts * Standarized.factors["SD"]) + Standarized.factors["Mean"]
high.limit.mts_no_s <- (high.limit.mts * Standarized.factors["SD"]) + Standarized.factors["Mean"]

# - remove log transformation from mortality
mt <- exp(1)^mt_no_s - 1/20000
mts <- exp(1)^mts_no_s - 1/20000

## - For the CI's we have again to remove log transformation on the all CI calculation, 
## - not individually on Ct, Qt and Cts
low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000

# - plot with the real values
min <- min(low.limit.mt, na.rm=T)
max <- max(high.limit.mt, na.rm=T)

#mt - darkgreen
plot(mt~results$mt$date, type="l", xlab="date", ylab="mortality", 
     main="Filtered mean - Country", ylim=c(min,max), lwd=2, col='darkgreen') 
lines(low.limit.mt~results$mt$date, col='darkgreen', lty = "dashed")
lines(high.limit.mt~results$mt$date, col='darkgreen', lty = "dashed") 
legend(x = "bottomright",                    
       legend = c("mt", "95% CI"),  
       col = c("darkgreen", "darkgreen"),   
       lwd = 1,                             
       lty = c(1,2),
       cex=0.7)


#mts (smoother) - blue
plot(mts~results$mts$date, type="l", xlab="date", ylab="mortality", 
     main="Smoothed mean - Country", ylim=c(min,max), lwd=2, col='blue') 
lines(low.limit.mts~results$mts$date, col='blue', lty = "dashed")
lines(high.limit.mts~results$mts$date, col='blue', lty = "dashed") 
legend(x = "bottomright",                    
       legend = c("mts", "95% CI"),  
       col = c("blue", "blue"),   
       lwd = 1,                             
       lty = c(1,2),
       cex=0.7)



## Region ----

regions <- c("A", "B", "C", "D", "E")

for (region in regions){
  
  # - Plot with log transformed and standardized values
  min <- min(results$mt[,region]-1.96*sqrt(results$Ct[,region]), na.rm=T)
  max <- max(results$mt[,region]+1.96*sqrt(results$Ct[,region]), na.rm=T)
  
  #mt - darkgreen
  mt=results$mt[,region]
  low.limit.mt <- (results$mt[,region])-1.96*sqrt(results$Ct[,region])
  high.limit.mt <- (results$mt[,region])+1.96*sqrt(results$Ct[,region])

  plot(mt~results$mt$date, type="l", xlab="date", ylab="log.mortality", 
       main=paste("Filtered mean - region: - ", region, sep=""), ylim=c(min,max), lwd=2, col='darkgreen') 
  lines(low.limit.mt~results$mt$date, col='darkgreen', lty = "dashed")
  lines(high.limit.mt~results$mt$date, col='darkgreen', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("mt", "95% CI"),  
         col = c("darkgreen", "darkgreen"),   
         lwd = 1,                             
         lty = c(1,2),
         cex=0.7)
  
  
  #mts (smoother) - blue
  mts=results$mts[,region]
  low.limit.mts <- (results$mts[,region])-1.96*sqrt(results$Cts[,region])
  high.limit.mts <- (results$mts[,region])+1.96*sqrt(results$Cts[,region])
  
  plot(mts~results$mts$date, type="l", xlab="date", ylab="log.mortality", 
       main=paste("Smoothed mean - region: - ", region, sep=""), ylim=c(min,max), lwd=2, col='blue') 
  lines(low.limit.mts~results$mts$date, col='blue', lty = "dashed")
  lines(high.limit.mts~results$mts$date, col='blue', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("mts", "95% CI"),  
         col = c("blue", "blue"),   
         lwd = 1,                             
         lty = c(1,2),
         cex=0.7)
  
  
  # - remove standardization from mortality
  mt_no_s <- (results$mt[,region] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  mts_no_s <- (results$mts[,region] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  ## - For the CI's we have to remove the standardization on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt_no_s <- (unlist(low.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mt_no_s <- (unlist(high.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  low.limit.mts_no_s <- (unlist(low.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mts_no_s <- (unlist(high.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  # - remove log transformation from mortality
  mt <- exp(1)^mt_no_s - 1/20000
  mts <- exp(1)^mts_no_s - 1/20000

  ## - For the CI's we have again to remove log transformation on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
  high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
  low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
  high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000
  
  # - plot with the real values
  min <- min(low.limit.mt, na.rm=T)
  max <- max(high.limit.mt, na.rm=T)
  
  #mt - darkgreen
  plot(mt~results$mt$date, type="l", xlab="date", ylab="mortality", 
       main=paste("Filtered mean - region: - ", region, sep=""), ylim=c(min,max), lwd=2, col='darkgreen') 
  lines(low.limit.mt~results$mt$date, col='darkgreen', lty = "dashed")
  lines(high.limit.mt~results$mt$date, col='darkgreen', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("mt", "95% CI"),  
         col = c("darkgreen", "darkgreen"),   
         lwd = 1,                             
         lty = c(1,2),
         cex=0.7)

  #mts - blue
  plot(mts~results$mts$date, type="l", xlab="date", ylab="mortality", 
       main=paste("Smoothed mean - region: - ", region, sep=""), ylim=c(min,max), lwd=2, col='blue') 
  lines(low.limit.mts~results$mts$date, col='blue', lty = "dashed")
  lines(high.limit.mts~results$mts$date, col='blue', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("mts", "95% CI"),  
         col = c("blue", "blue"),   
         lwd = 1,                             
         lty = c(1,2),
         cex=0.7)
}


## Farms ----

farms <- colnames(D[2:ncol(D)])
#farm="451_A"
#farm="456_C"

for (farm in farms){
  
  # - Plot with log transformed and standardized values
  min <- min(results$ft[,farm]-1.96*sqrt(results$Qt[,farm]), na.rm=T)
  max <- max(results$ft[,farm]+1.96*sqrt(results$Qt[,farm]), na.rm=T)
  
  #mt - darkgreen
  mt=results$mt[,farm]
  low.limit.mt <- (results$mt[,farm])-1.96*sqrt(results$Ct[,farm])
  high.limit.mt <- (results$mt[,farm])+1.96*sqrt(results$Ct[,farm])

  plot(mt~results$mt$date, type="l", xlab="date", ylab="log.mortality", 
       main=paste("Filtered mean - farm: - ", farm, sep=""), ylim=c(min,max), lwd=2, col='darkgreen') 
  lines(D[,farm]~results$mt$date, lwd=2)
  lines(low.limit.mt~results$mt$date, col='darkgreen', lty = "dashed")
  lines(high.limit.mt~results$mt$date, col='darkgreen', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("mt", "Obs", "95% CI"),  
         col = c("darkgreen", "black", "darkgreen"),   
         lwd = 1,                             
         lty = c(1,1,2),
         cex=0.7)
  
  #mts (smoother) - blue
  mts=results$mts[,farm]
  low.limit.mts <- (results$mts[,farm])-1.96*sqrt(results$Cts[,farm])
  high.limit.mts <- (results$mts[,farm])+1.96*sqrt(results$Cts[,farm])
  
  plot(mts~results$mts$date, type="l", xlab="date", ylab="log.mortality", 
       main=paste("Smoothed mean - farm: - ", farm, sep=""), ylim=c(min,max), lwd=2, col='blue') 
  lines(D[,farm]~results$mts$date, lwd=2)
  lines(low.limit.mts~results$mts$date, col='blue', lty = "dashed")
  lines(high.limit.mts~results$mts$date, col='blue', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("mts", "Obs", "95% CI"),  
         col = c("blue", "black", "blue"),   
         lwd = 1,                             
         lty = c(1,1,2),
         cex=0.7)

  #ft - red
  ft=results$ft[,farm]
  low.limit.ft <- (results$ft[,farm])-1.96*sqrt(results$Qt[,farm])
  high.limit.ft <- (results$ft[,farm])+1.96*sqrt(results$Qt[,farm])

  plot(ft~results$ft$date, type="l", xlab="date", ylab="log.mortality", 
       main=paste("Forecasts - farm: - ", farm, sep=""), ylim=c(min,max), lwd=2, col='red') 
  lines(D[,farm]~results$ft$date, lwd=2)
  lines(low.limit.ft~results$ft$date, col='red', lty = "dashed")
  lines(high.limit.ft~results$ft$date, col='red', lty = "dashed") 
  legend(x = "bottomright",                    
         legend = c("ft", "Obs", "95% CI"),  
         col = c("red", "black", "red"),   
         lwd = 1,                             
         lty = c(1,1,2),
         cex=0.7)
  
  
  # - remove standardization from mortality
  mt_no_s <- (results$mt[,farm] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  ft_no_s <- (results$ft[,farm] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  mts_no_s <- (results$mts[,farm] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  ## - For the CI's we have to remove the standardization on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt_no_s <- (unlist(low.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mt_no_s <- (unlist(high.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  low.limit.ft_no_s <- (unlist(low.limit.ft) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.ft_no_s <- (unlist(high.limit.ft) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  low.limit.mts_no_s <- (unlist(low.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mts_no_s <- (unlist(high.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  # - remove log transformation from mortality
  Obs <- exp(1)^(Test.set[,farm]) - 1/20000
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
  
  #mt - darkgreen
  plot(mt~results$mt$date, type="l", xlab="date", ylab="mortality", 
       main=paste("Filtered mean - farm: - ", farm, sep=""), ylim=c(min,max), lwd=2, col='darkgreen') 
  lines(Obs~results$mt$date, lwd=2)
  lines(low.limit.mt~results$mt$date, col='darkgreen', lty = "dashed")
  lines(high.limit.mt~results$mt$date, col='darkgreen', lty = "dashed") 
  legend(x = "topright",                    
         legend = c("mt", "Obs", "95% CI"),  
         col = c("darkgreen", "black", "darkgreen"),   
         lwd = 1,                             
         lty = c(1,1,2),
         cex=0.7)
  
  
  #mts (smoother) - blue
  plot(mts~results$mts$date, type="l", xlab="date", ylab="mortality", 
       main=paste("Smoothed mean - farm: - ", farm, sep=""), ylim=c(min,max), lwd=2, col='blue') 
  lines(Obs~results$mts$date, lwd=2)
  lines(low.limit.mts~results$mts$date, col='blue', lty = "dashed")
  lines(high.limit.mts~results$mts$date, col='blue', lty = "dashed") 
  legend(x = "topright",                    
         legend = c("mts", "Obs", "95% CI"),  
         col = c("blue", "black", "blue"),   
         lwd = 1,                             
         lty = c(1,1,2),
         cex=0.7)
  
  #ft - red
  plot(ft~results$ft$date, type="l", xlab="date", ylab="mortality", 
       main=paste("Forecasts - farm: - ", farm, sep=""), ylim=c(min,max), lwd=2, col='red') 
  lines(Obs~results$ft$date, lwd=2)
  lines(low.limit.ft~results$mts$date, col='red', lty = "dashed")
  lines(high.limit.ft~results$mts$date, col='red', lty = "dashed") 
  legend(x = "topright",                    
         legend = c("ft", "Obs", "95% CI"),  
         col = c("red", "black", "red"),   
         lwd = 1,                             
         lty = c(1,1,2),
         cex=0.7)
}


