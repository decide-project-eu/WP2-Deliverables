#### Functions for Multi-level DLM - DECIDE Deliverable 2.2


# Function to estimate the mu0 (initial parameter vector) ----
# - D is the data set, containing the time series we want to model
# - D.months contains the information about the start and duration of the production cycles
# - expected.start.time is the observation time when we expect the model to start (usually 0 or 1)
# - regions is a vector with the name of the regions
get.m0 <- function(D, D.months, expected.start.time, regions){
  
  ### Elements for the initial parameter vector (mu0)
  mu0 <- c()
  Spline.list <- list()
  
  # Remove outliers, based on overall mean and moving SD
  all.values <- unlist(D[,c(2:ncol(D))])
  ylim <- range(na.omit(all.values))
  Mean <- mean(na.omit(all.values))
  SD <- sd(na.omit(all.values))
  upper <- Mean + 3*SD #99.7% of data occurs within 3 standard deviations of the mean within a normal distribution
  lower <- Mean - 3*SD
  for(farm in colnames(D[2:ncol(D)])){
    remove.i <- which(D[,farm] > upper | D[,farm] < lower)
    D[remove.i, farm] <- NA
  }

  #Country
  m0.C <- mean(unlist(D[,c(2:ncol(D))]), na.rm=T)
  mu0 <- rbind(mu0, m0.C)
  row.names(mu0) <- "Country"
  
  #Regions
  regions.means <- c()
  
  for (region in regions){ 
    farms.per.region <- grep(region, names(D), value=TRUE)
    m0.R <- mean(unlist(D[,farms.per.region]), na.rm=T)
    regions.means <- c(regions.means, m0.R)
    names(regions.means)[length(regions.means)] <- region
    m0.R.C <- m0.R - m0.C
    mu0 <- rbind(mu0, m0.R.C)
    row.names(mu0)[nrow(mu0)] <- region
  }
  
  #Farms
  for(farm in colnames(D[2:ncol(D)])){ 
    
    mu <- c()
    y.i <- D[,farm]
    y <- y.i[-which(is.na(y.i))]
    
    x.i <- D.months[,farm]
    x <- x.i[-which(is.na(y.i))] # because what is NA in mortality (not known or fallow period) those months have to be removed
    
    Spline = smooth.spline(x = x, y = y)
    
    plot(y~x, xlab="month", ylab=farm, main = paste("Spline", '=', farm))
    lines(Spline, col='red', lwd=3)
    
    # - Save the spline function - we will need it for the Gt matrix
    spline.name <- paste(farm, '_Spline', sep='')
    Spline.list[[length(Spline.list)+1]] <- Spline
    names(Spline.list)[length(Spline.list)] <- spline.name
    
    # - Get the region mean (of the corresponding farm)
    region <- sapply(strsplit(farm, '_'), `[`, 2) 
    m0.R <- regions.means[region]
    
    # - Finalize and save mu
    pred0 <- predict(object = Spline, x = (expected.start.time-1))$y
    m0.F.R.C <- pred0 - m0.R - m0.C
    mu <- c(m0.F.R.C, 1)
    
    #give names to columns and rows
    mu.names <- c(farm, paste('d.', farm, sep=''))      
    mu <- matrix(mu)
    row.names(mu) <- mu.names
    
    # - Initial parameter vector for all variables (mu0)
    mu0 <- rbind(mu0, mu)
    colnames(mu0) <- "mu0"
  }
  return(list('mu0'=mu0, 'Spline.list'=Spline.list))
}



# Function to estimate C0 (prior variance) ----
# - D is the data set, containing the time series we want to model
# - D.months contains the information about the start and duration of the production cycles
# - expected.start.time is the observation time when we expect the model to start (usually 0 or 1)
# - regions is a vector with the name of the regions
get.C0 <- function(D, D.months, expected.start.time, regions){
  
  # Remove outliers, based on overall mean and moving SD
  all.values <- unlist(D[,c(2:ncol(D))])
  ylim <- range(na.omit(all.values))
  Mean <- mean(na.omit(all.values))
  SD <- sd(na.omit(all.values))
  upper <- Mean + 3*SD #99.7% of data occurs within 3 standard deviations of the mean within a normal distribution
  lower <- Mean - 3*SD
  for(farm in colnames(D[2:ncol(D)])){
    remove.i <- which(D[,farm] > upper | D[,farm] < lower)
    D[remove.i, farm] <- NA
  }
  
  # Make the empty C0 matrix
  C.all <- list()
  n.rows.cols <- 1+5+length(colnames(D[2:ncol(D)]))*2
  C0 <- as.data.frame(matrix(0, ncol = n.rows.cols, nrow = n.rows.cols))
  
  # Give row and column names to C0
  country.regions.names <- c("Country", regions)
  level.trend.names <- c()
  for(farm in colnames(D[2:ncol(D)])){
    level.trend.names <- c(level.trend.names, farm, paste('d.', farm, sep=''))      
  }
  row.col.names <- c(country.regions.names, level.trend.names) 
  row.names(C0) <- row.col.names
  colnames(C0) <- row.col.names
  
  # Make an empty table for saving country, regions and farms means
  means <- as.data.frame(matrix(NA, ncol=1+length(regions)+length(colnames(D[2:ncol(D)])), nrow=length((expected.start.time):(expected.start.time+6))))
  row.names(means) <- (expected.start.time):(expected.start.time+6)
  colnames(means) <- c("Country", regions, colnames(D[2:ncol(D)]))
  
  
  # Elements for the prior variance (C0)
  #Country
  
  # - Create table with the relevant information
  D.A <- data.frame(Country=as.vector(unlist(D[,c(2:ncol(D))])), months=as.vector(unlist(D.months[,c(2:ncol(D.months))])))
  #only consider the first 6 months
  D.A <- subset(D.A, D.A[,"months"] %in% (expected.start.time):(expected.start.time+6))
  
  # - Get the mean of the country for the first 7 months
  for(i in (expected.start.time):(expected.start.time+6)){
    D.A.i <- subset(D.A, D.A$months == i)
    means[row.names(means)[row.names(means)==i], "Country"] <- mean(D.A.i$Country, na.rm=T)
  }

  # - Finalize and save C0
  C <- var(means$Country)
  C0["Country", "Country"] <- C
  
  # - Save all C matrices to see where there is NA values and replace them afterwards
  C.all[["Country"]] <- C
  
  
  #Regions
  for (region in regions){
    
    # - Get the farms that belong to the current region
    indx.region <- grepl(region, colnames(D[,c(2:ncol(D))]))
    farms.region <- colnames(D[,c(2:ncol(D))])[indx.region]
    
    # - Create table with the relevant information
    D.A <- data.frame(region=as.vector(unlist(D[,farms.region])), months=as.vector(unlist(D.months[,farms.region])))
    #only consider the first 7 months
    D.A <- subset(D.A, D.A[,"months"] %in% (expected.start.time):(expected.start.time+6))
    colnames(D.A)[1] <- region
    
    # - Get the mean of each region for the first 7 months
    for(i in (expected.start.time):(expected.start.time+6)){
      D.A.i <- subset(D.A, D.A$months == i)
      means[row.names(means)[row.names(means)==i], region] <- mean(D.A.i[, region], na.rm=T)
    }
    
    # - Get the difference between the means of the region and the means of the country
    R.C <- means[, region] - means[, "Country"]
    
    # - Finalize and save C0
    C <- var(R.C)
    C0[region, region] <- C
    
    # - Save all C matrices to see where there is NA values and replace them afterwards
    C.all[[region]] <- C
  }
  
  
  #Farms
  for(farm in colnames(D[2:ncol(D)])){
    
    # - Create table with the relevant information
    D.A <- data.frame(Farm=D[,farm], months=D.months[,farm])
    #only consider the first 6 months
    D.A <- subset(D.A, D.A[,"months"] %in% (expected.start.time):(expected.start.time+6))
    colnames(D.A)[1] <- farm
    
    # - Get the mean of each farm for the first 7 months
    for(i in (expected.start.time):(expected.start.time+6)){
      D.A.i <- subset(D.A, D.A$months == i)
      means[row.names(means)[row.names(means)==i], farm] <- mean(D.A.i[, farm], na.rm=T)
    }
    
    # - Get the region which this farm belongs to
    region <- sapply(strsplit(farm, '_'), `[`, 2)
    
    # - Get the difference between the means of the farm, region and the country
    F.R.C <- means[, farm] - means[, region] - means[, "Country"]
    
    # - Difference between the present observation and the previous observation (deviations)
    Diff <- diff(F.R.C)
    F.R.C.Diff <- c(NA, Diff)
    
    # - Combine the differences and the deviations
    df <- data.frame(F.R.C=F.R.C, F.R.C.Diff=F.R.C.Diff)
    df <- na.omit(df)
    
    # - Finalize and save C0
    C <- cov(df[,c("F.R.C", "F.R.C.Diff")])
    colnames(C) <- c(farm, paste("d.", farm, sep=""))
    row.names(C) <- c(farm, paste("d.", farm, sep=""))
    C0[c(farm, paste("d.", farm, sep="")), c(farm, paste("d.", farm, sep=""))] <- C
    C0 <- as.matrix(C0)

    # - Save all C matrices to see where there is NA values and replace them afterwards
    C.all[[farm]] <- C
  }
  return(list('C0'=C0, 'C.all'=C.all))
}

  

# Function to estimate the Ft (design matrix) ----
# - D is the data set, containing the time series we want to model
# - regions is a vector with the name of the regions
get.Ft <- function(D, regions){
  
  # Make the empty Ft matrix
  n.rows <- 1+5+length(colnames(D[2:ncol(D)]))*2
  n.cols <- length(colnames(D[2:ncol(D)]))
  Ft <- as.data.frame(matrix(0, ncol = n.cols, nrow = n.rows))
  
  # Give row and column names to Ft
  # - Rows
  country.regions.names <- c("Country", regions)
  level.trend.names <- c()
  for(farm in colnames(D[2:ncol(D)])){
    level.trend.names <- c(level.trend.names, farm, paste('d.', farm, sep=''))      
  }
  row.names <- c(country.regions.names, level.trend.names) 
  row.names(Ft) <- row.names
  # - Columns
  colnames(Ft) <- colnames(D[2:ncol(D)])
  
  # Ft matrix
  for(farm in colnames(Ft)){ 
    Ft["Country", farm] <- 1
    region <- sapply(strsplit(farm, '_'), `[`, 2) 
    Ft[region, farm] <- 1
    Ft[farm, farm] <- 1
  }
  return(as.matrix(Ft))
}
  


# Function to estimate the Gt (system matrix) ----
# - D is the data set, containing the time series we want to model
# - D.months contains the information about the start and duration of the production cycles
# - i is the time step you are modelling in the DLM, you can also give it equal to NA if you want to use dates instead
# - Date is the date you are currently modelling in the DLM
# - Date_1 is the previously date to what you are currently modelling in the DLM
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector with the name of the regions
get.Gt <- function(D, D.months, i, Date, Date_1, Spline.list, regions){
  
  # Make the empty Gt matrix
  n.rows.cols <- 1+5+length(colnames(D[2:ncol(D)]))*2
  Gt <- as.data.frame(diag(1, ncol = n.rows.cols, nrow = n.rows.cols))
  
  # Give row and column names to Gt
  country.regions.names <- c("Country", regions)
  level.trend.names <- c()
  for(farm in colnames(D[2:ncol(D)])){ 
    level.trend.names <- c(level.trend.names, farm, paste('d.', farm, sep=''))      
  }
  row.col.names <- c(country.regions.names, level.trend.names) 
  row.names(Gt) <- row.col.names
  colnames(Gt) <- row.col.names
  
  
  # Elements for the system matrix (Gt)
  
  for(farm in colnames(D[2:ncol(D)])){ 
    
    # - get the relevant spline function
    spline.name <- paste(farm, 'Spline', sep='_') #1 spline per farm
    Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
    
    # - get the date index
    i <- which(unique(D$date) == Date)
    i_1 <- which(unique(D$date) == Date_1)
    
    if(isTRUE(length(i)!=0)){
      
      time <- D.months[i,farm] #The observation time of the current observations
      time_1 <- D.months[i_1,farm] #The observation time of the previous observations
      if(length(time_1) == 0){
        time_1 <- 0
      }
      
      # - get trend for Gt matrix
      if(!is.na(time)){
        pred <- predict(object = Spline, x = time)$y
        if(isTRUE(time != time_1)){
          pred_1 <- predict(object = Spline, x = time_1)$y
          Trend <- pred - pred_1
        }else{
          Trend <- pred 
        }
      }else{
        Trend <- 0
      }
      
      Gt[farm, paste("d.", farm, sep="")] <- Trend
      
    }else{
      Gt[farm, paste("d.", farm, sep="")] <- 0
    }
  }
  return(as.matrix(Gt))
}
  
  

# Function to estimate the V (observation variance) ----
# - D is the data set, containing the time series we want to model
# - countryVar is the variance component of country
# - regionVar is the variance component of region
# - farmVar is the variance component of farm
get.Vt <- function(D, countryVar, regionVar, farmVar){
  
  # Make the empty Vt matrix
  n.rows.cols <- length(colnames(D[2:ncol(D)]))
  Vt <- as.data.frame(matrix(0, ncol = n.rows.cols, nrow = n.rows.cols))
  
  # Give row and column names to Ft
  rownames(Vt) <- colnames(D[2:ncol(D)])
  colnames(Vt) <- colnames(D[2:ncol(D)])
  
  for(farm in rownames(Vt)){
    
    #Country
    Vt[farm, ] <- countryVar
    
    #Region
    region <- sapply(strsplit(farm, '_'), `[`, 2)
    indx.region <- grepl(region, rownames(Vt))
    farms.region <- rownames(Vt)[indx.region]
    Vt[farm, farms.region] <- countryVar + regionVar
    
    #Farm
    Vt[farm, farm] <- countryVar + regionVar + farmVar
  }
  Vt <- as.matrix(Vt)
  return(Vt)
}



# Function to initialized mt ----
# - mt is the filtered mean calculated my the DLM
# - mu0 is the initial parameter vector
# - D.months contains the information about the start and duration of the production cycles
# - Date is the date you are currently modelling in the DLM
## when a new production cycle starts the mt for that farm should be initialized (=mu0) 
get.mt.hierarchical <- function(mt, mu0, D.months, Date){
  
  #Observation vector only for that Date
  Date.set.months <- D.months[which(unique(D.months$date) == Date),2:ncol(D.months)]
  
  #Farms that started a new production cycle on that Date
  Date.first.set <- which(Date.set.months==0)
  farms.first <- colnames(Date.set.months[Date.first.set])
  
  mt.h <- mt
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){ 
      
      #get level and trends names for those farms
      farm.name.mu0.l <- farm
      farm.name.mu0.t <- paste("d", farm, sep=".")
      
      #replace mt for those farms by mu0 (initialize)
      mt.h[c(farm.name.mu0.l, farm.name.mu0.t),] <- as.matrix(mu0[c(farm.name.mu0.l, farm.name.mu0.t),])
    }
  }
  return(mt.h)
}



# Function to initialized Ct ----
# - Ct is the filtered variance calculated my the DLM
# - C0 is the prior variance
# - D.months contains the information about the start and duration of the production cycles
# - Date is the date you are currently modelling in the DLM
## when a new production cycle starts the Ct for that farm should be initialized (=C0) and
## covariances should have the same correlation but different magnitudes
get.Ct.hierarchical <- function(Ct, C0, D.months, Date){
  
  #Observation vector only for that Date
  Date.set.months <- D.months[which(unique(D.months$date) == Date),2:ncol(D.months)]
  
  #Farms that started a new production cycle on that Date
  Date.first.set <- which(Date.set.months==0)
  farms.first <- colnames(Date.set.months[Date.first.set])
  
  Ct.h <- Ct
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){
      
      #get level and trends names for those farms
      farm.name.C0.l <- farm
      farm.name.C0.t <- paste("d", farm, sep=".")
      
      #replace Ct for those farms by C0 (initialize)
      # - replace 4 "diagonal" values
      Ct.h[c(farm.name.C0.l, farm.name.C0.t),c(farm.name.C0.l, farm.name.C0.t)] <- 
        as.matrix(C0[c(farm.name.C0.l, farm.name.C0.t),c(farm.name.C0.l, farm.name.C0.t)])
      
      #replace covariances
      # - get variance of C0 and Ct
      C0.a <- as.matrix(C0[farm.name.C0.l, farm.name.C0.l]) #variance of level
      C0.b <- as.matrix(C0[farm.name.C0.t, farm.name.C0.t]) #variance of trend
      Ct.a <- as.matrix(Ct[farm.name.C0.l, farm.name.C0.l]) #variance of level
      Ct.b <- as.matrix(Ct[farm.name.C0.t, farm.name.C0.t]) #variance of trend
      
      # - start replacing the rows and the columns that have the level (a)
      # - first rows and columns upper left side (1)
      cov.a1 <- t(as.matrix(Ct.h[farm.name.C0.l, 1:((which(colnames(Ct.h)==farm.name.C0.l))-1)]))
      rownames(cov.a1) <- farm.name.C0.l
      
      for(i in 1:length(cov.a1)){
        Ct.h[farm.name.C0.l,i] <- (cov.a1[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
        Ct.h[i,farm.name.C0.l] <- (cov.a1[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
      }
      
      # - second rows and columns down right side (2)
      if(farm.name.C0.l != rownames(C0)[(dim(C0)[1])-1]){
        cov.a2 <- t(as.matrix(Ct.h[farm.name.C0.l, ((which(colnames(Ct.h)==farm.name.C0.t))+1):dim(Ct.h)[1]]))
        rownames(cov.a2) <- farm.name.C0.l
        
        for(i in 1:length(cov.a2)){
          row.name <- colnames(cov.a2)[i] 
          Ct.h[farm.name.C0.l,row.name] <- (cov.a2[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
          Ct.h[row.name,farm.name.C0.l] <- (cov.a2[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
        }
      }
      
      # - now replace the rows and the columns that have the trend (b)
      # - first rows and columns upper left side (1)
      cov.b1 <- t(as.matrix(Ct.h[farm.name.C0.t, 1:((which(colnames(Ct.h)==farm.name.C0.l))-1)]))
      rownames(cov.b1) <- farm.name.C0.t
      
      for(i in 1:length(cov.b1)){
        Ct.h[farm.name.C0.t,i] <- (cov.b1[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
        Ct.h[i,farm.name.C0.t] <- (cov.b1[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
      }
      
      # - second rows and columns down right side (2)
      if(farm.name.C0.l != rownames(C0)[(dim(C0)[1])-1]){ 
        cov.b2 <- t(as.matrix(Ct.h[farm.name.C0.t, ((which(colnames(Ct.h)==farm.name.C0.t))+1):dim(Ct.h)[1]]))
        rownames(cov.b2) <- farm.name.C0.t
        
        for(i in 1:length(cov.b2)){
          row.name <- colnames(cov.b2)[i] 
          Ct.h[farm.name.C0.t,row.name] <- (cov.b2[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
          Ct.h[row.name,farm.name.C0.t] <- (cov.b2[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
        }
      }
    }
  }
  return(Ct.h)
}



# Function of the filtering algorithm (DLM) with discount factors ----
# - D is the data set, containing the time series we want to model
# - D.months contains the information about the start and duration of the production cycles
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - countryVar is the variance component of country
# - regionVar is the variance component of region
# - farmVar is the variance component of farm
# - deltas (country, region and farm - in this order) are the discount factors used to create system variance-covariance matrix (W)
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector with the name of the regions
runDiscountDLM <- function(D, D.months, mu0, C0, countryVar, regionVar, farmVar, deltas, Spline.list, regions){
  
  n <- 1:length(unique(D$date))
  first.date <- as.Date(D$date)[1]
  
  Yt.list <- list()
  at.list <- list()		 
  Rt.list <- list()
  ft.list <- list()
  Qt.list <- list()
  At.list <- list()
  et.list <- list()
  ut.list <- list()
  mt.list <- list()
  Ct.list <- list()
  Ft.list <- list()
  VSE.list <- list()
  Gt.list <- list()
  fullFt.list <- list()
  
  # Create list for discount groups
  discount.groups <- list()
  discount.groups[[1]] <- 1             #country
  discount.groups[[2]] <- 2:6           #5 regions
  discount.groups[[3]] <- 7:dim(C0)[1]  #all farms (levels and trends)
  
  mt <- mu0				
  Ct <- C0
  
  # Make sure Ct is symmetrical
  Ct <- (Ct + t(Ct))/2
  
  for(i in n){
    
    Date <- sort(unique(D$date))[i]
    Date_1 <- sort(unique(D$date))[i-1]
    if(length(Date_1)==0){
      Date_1 <- 0
    }
    
    # Get the observation vector
    Yt <- D[which(unique(D$date) == Date),2:ncol(D)]
    
    # Get the observational variance (Vt)
    Vt <- get.Vt(D, countryVar, regionVar, farmVar)
    
    # Get Ft given the current Yt 
    Ft <- get.Ft(D, regions)
    
    # Get Gt - independent of the current Yt
    Gt <- get.Gt(D, D.months, i=NA, Date, Date_1, Spline.list, regions)
    
    
    # Remove missing values from Yt, Vt and Ft
    missing = which(is.na(Yt))
    fullYt <- Yt
    fullFt <- Ft
    fullVt <- Vt
    
    # - If length of missing is > 0 there is at least one missing
    if (length(missing) > 0) {
      
      # remove from Yt
      a <- colnames(fullYt)[-missing]
      Yt <- t(as.matrix(fullYt[-missing]))
      
      # Remove from Ft
      Ft <- as.matrix(fullFt[, -missing])
      
      # Remove from Vt
      Vt <- fullVt[-missing, -missing]
    }
    
    # Start with the Kalman filter
    at <- Gt %*% mt                        
    Pt <- Gt %*% Ct %*% t(Gt)              
    
    # Define system variance
    Wt <- matrix(0, ncol = length(C0[,1]), nrow = length(C0[,1]))
    for (b in 1:3) {
      Wt[discount.groups[[b]], discount.groups[[b]]] = ((1-deltas[b])/deltas[b]) * Pt[discount.groups[[b]], discount.groups[[b]]]
    }
    
    # If all observations are NA
    if (length(Yt) == 0) {
      
      ## Initialize mt and Ct when a new production cycle starts or when starts to have information for one farm
      if(Date > first.date){
        mt <- get.mt.hierarchical(mt, mu0, D.months, Date)
        Ct <- get.Ct.hierarchical(Ct, C0, D.months, Date)
      }
      
      ## Return to Kalman Filter
      Rt = Pt + Wt                                   # Prior Variance
      ft = t(fullFt) %*% at                          # One-step Forecast mean
      Qt = t(fullFt) %*% Rt %*% fullFt + fullVt      # One-step Forecast variance
      At = NA                                        # Adaptative Coef. matrix
      et = NA	                                       # One-step forecast error
      ut = NA                                        # Standardized forecast error
      
      # - update the parameter vector and variance matrix
      mt <- at                                       # Filtered mean
      Ct <- Rt	                                     # Filtered variance
      
    }else{ #if we observe at least one farm
      
      ## Initialize mt and Ct when a new production cycle starts or when starts to have information for one farm
      if(Date > first.date){
        mt <- get.mt.hierarchical(mt, mu0, D.months, Date)
        Ct <- get.Ct.hierarchical(Ct, C0, D.months, Date)
      }
      
      ## Return to Kalman Filter
      Rt = Pt + Wt                                   # Prior Variance
      ft = t(Ft) %*% at                              # One-step Forecast mean
      ft2 = t(fullFt) %*% at                         # Save this one to have estimates of ft when we don't have observations
      Qt = t(Ft) %*% Rt %*% Ft + Vt                  # One-step Forecast variance
      Qt2 = t(fullFt) %*% Rt %*% fullFt + fullVt     # Save this one to have estimates of Qt when we don't have observations
      At = Rt %*% Ft %*% solve(Qt)                   # Adaptative Coef. matrix
      et <- Yt  - ft	                               # One-step forecast error
      ut <- et / sqrt(diag(Qt))                      # Standardized forecast error
      
      # - update the parameter vector and variance matrix
      mt <- at + At %*% et                           # Filtered mean
      Ct <- Rt - At  %*% Qt %*% t(At)	               # Filtered variance
    }
    
    # Make sure Ct is symmetrical
    Ct <- (Ct + t(Ct))/2
    
    # Save the values in lists
    Yt.list[[i]] <- Yt
    at.list[[i]] <- at
    Rt.list[[i]] <- Rt
    ft.list[[i]] <- ft2
    Qt.list[[i]] <- Qt2
    At.list[[i]] <- At
    et.list[[i]] <- et
    ut.list[[i]] <- ut
    mt.list[[i]] <- mt
    Ct.list[[i]] <- Ct
    Ft.list[[i]] <- t(Ft)
    Gt.list[[i]] <- Gt
    fullFt.list[[i]] <- t(fullFt)
  }
  
  return(list(
    Yt=Yt.list,
    at=at.list,
    Rt=Rt.list,
    ft=ft.list,
    Qt=Qt.list,
    At=At.list,
    et=et.list,
    ut=ut.list,
    mt=mt.list,
    Ct=Ct.list,
    F=Ft.list,
    Gt.list=Gt.list,
    fullFt.list=fullFt.list
  ))
}



# Function of the Kalman Smoother ----
# - res is a result returned from the filter (DLM)
runSmoother <- function(res) {
  
  n = length(res$mt)
  p = length(res$mt[[1]])
  mts <- array(NA,dim=c(p,1,n));
  Cts <- array(NA,dim=c(p,p,n));
  
  # Put last value equal to filtered
  mts[,,n] <- res$mt[[n]]
  Cts[,,n] <- res$Ct[[n]]
  
  # These are useful
  Bt <- array(NA,dim=c(p,p,n))
  Lt <- array(NA,dim=c(p,p,n));  
  
  # Iterate backwards over days
  for(i in ((n-1):1))   {
    
    # Get Gt
    Gt <- res$Gt.list[[i+1]]
    
    res$R[[i+1]] <- as.matrix(res$R[[i+1]])
    
    Bt[,,i] <- as.matrix( res$Ct[[i]] %*% t(Gt) %*% solve(res$R[[i+1]]) )
    mts[,,i] <- res$mt[[i]] + Bt[,,i] %*% (mts[,,i+1] - res$a[[i+1]])
    Cts[,,i] <- as.matrix( res$C[[i]] + Bt[,,i] %*% (Cts[,,i+1] - res$R[[i+1]]) %*% t(Bt[,,i]) )
    
    mts[,,i] <- as.matrix(mts[,,i])
  }
  # give names
  rownames(mts) <- rownames(res$mt[[i]])
  colnames(mts) <- "mts"
  rownames(Cts) <- rownames(res$Ct[[i]])
  colnames(Cts) <- colnames(res$Ct[[i]])
  
  # Now when we are at it: Find L and store it for the EM algorithm
  for(i in ((n):2))  {
    Lt[,,i] <- Cts[,,i] + Gt%*%Cts[,,i-1]%*%t(Gt) - Cts[,,i]%*%t(Bt[,,i-1]) - Bt[,,i-1]%*%Cts[,,i]
  }
  rownames(Lt) <- rownames(res$Ct[[i]])
  colnames(Lt) <- colnames(res$Ct[[i]])
  
  return(list(mts=mts,
              Cts=Cts,
              Lt=Lt,
              D=res$D));
}



## Function to extract the relevant information ----
# - res is a result returned from the filter (DLM)
# - smot is a result returned from the smoother
# - D is the data set, containing the time series we want to model
# - H.w.country is in case you want to model the DLM with harmonic waves - here we don't so we put equal to NULL 
# - regions is a vector with the name of the regions
extract.res <- function(res, smot, D, H.w.country=NULL, regions){
  
  results.list <- list()
  
  # Get the filtered mean (mt)
  ## create empty data frame for mt results
  df.mt <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
  df.mt[,1] <- D$date
  colnames(df.mt) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## for country
  # - relevant vector to extract the country level
  vector <- rep(0, dim(res$mt[[1]])[1])
  names(vector) <- rownames(res$mt[[1]])
  
  if(is.null(H.w.country)){
    vector[1] <- 1
  }else{
    Fti <- 1
    for(n in 1:H.w.country){
      Fti <- c(Fti, c(1,0))
    }
    vector[c(1:(1+2*H.w.country))] <- Fti
  }
  vector.country <- vector #save vector for regions
  
  # - extract country level
  for(i in 1:length(res$mt)){ 
    c <- vector %*% res$mt[[i]]
    df.mt[i, "Country"] <- c
  }
  
  
  ## for region (sum of country and region)
  for(region in regions){
    
    # - relevant vector to extract the region level
    vector <- vector.country
    vector[region] <- 1
    
    # - extract region level
    for(i in 1:length(res$mt)){ 
      r <- vector %*% res$mt[[i]]
      df.mt[i, region] <- r
    }
  }
  
  
  ## for farms (sum of country, region and farm (use Ftfull))
  for(i in 1:length(res$mt)){ 
    f <- res$fullFt.list[[i]] %*% res$mt[[i]]
    df.mt[i, rownames(f)] <- t(f)
  }
  results.list[["mt"]] <- df.mt
  
  
  # Get the forecasts (ft)
  df.ft <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
  df.ft[,1] <- D$date
  colnames(df.ft) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  for(i in 1:length(res$ft)){ 
    for(a in rownames(res$ft[[i]])){
      df.ft[i, a] <- res$ft[[i]][a,1]
    }
  }
  results.list[["ft"]] <- df.ft
  
  
  # Get the raw forecasts errors (et)
  df.et <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
  df.et[,1] <- D$date
  colnames(df.et) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  
  for(i in 1:length(res$et)){ 
    for(a in rownames(res$et[[i]])){
      df.et[i, a] <- res$et[[i]][a,1]
    }
  }
  results.list[["et"]] <- df.et
  
  
  # Get the standardized forecasts errors (ut)
  df.ut <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
  df.ut[,1] <- D$date
  colnames(df.ut) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  for(i in 1:length(res$ut)){ 
    for(a in rownames(res$ut[[i]])){
      df.ut[i, a] <- res$ut[[i]][a,1]
    }
  }
  results.list[["ut"]] <- df.ut
  
  
  # Get the filtered variance (Ct) 
  ## create empty data frame for Ct results
  df.Ct <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
  df.Ct[,1] <- D$date
  colnames(df.Ct) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## for country
  # - relevant vector to extract the country level
  vector <- rep(0, dim(res$Ct[[1]])[1])
  names(vector) <- rownames(res$Ct[[1]])
  
  if(is.null(H.w.country)){
    vector[1] <- 1
  }else{
    Fti <- 1
    for(n in 1:H.w.country){
      Fti <- c(Fti, c(1,0))
    }
    vector[c(1:(1+2*H.w.country))] <- Fti
  }
  vector.country <- vector
  
  # - extract country level
  for(i in 1:length(res$Ct)){ 
    c <- vector %*% res$Ct[[i]] %*% vector
    df.Ct[i, "Country"] <- c
  }
  
  
  ## for region (considering country and region)
  for(region in regions){
    
    # - relevant vector to extract the region level
    vector <- vector.country
    vector[region] <- 1
    
    # - extract region level
    for(i in 1:length(res$Ct)){ 
      r <- vector %*% res$Ct[[i]] %*% vector 
      df.Ct[i, region] <- r
    }
  }
  
  
  ## for farms (considering country, region and farm (use Ftfull))
  for(i in 1:length(res$Ct)){ 
    f <- res$fullFt.list[[i]] %*% res$Ct[[i]] %*% t(res$fullFt.list[[i]])
    for(farm in colnames(f)){ 
      df.Ct[i, farm] <- f[farm, farm]
      
    }
  }
  results.list[["Ct"]] <- df.Ct
  
  
  # Get the forecast variance (Qt)
  df.Qt <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
  df.Qt[,1] <- D$date
  colnames(df.Qt) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  
  for(i in 1:length(res$Qt)){ 
    for(a in rownames(res$Qt[[i]])){  
      df.Qt[i, a] <- res$Qt[[i]][a,a]
    }
  }
  results.list[["Qt"]] <- df.Qt
  
  
  # Get the smoothed values
  if (!is.null(smot)) {
    
    # Get the smoothed mean (mts)
    ## create empty data frame for mts results
    df.mts <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
    df.mts[,1] <- D$date
    colnames(df.mts) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
    
    ## for country
    # - relevant vector to extract the country level
    vector <- rep(0, dim(smot$mts)[[1]])
    names(vector) <- names(smot$mts[,,1])
    
    if(is.null(H.w.country)){
      vector[1] <- 1
    }else{
      Fti <- 1
      for(n in 1:H.w.country){
        Fti <- c(Fti, c(1,0))
      }
      vector[c(1:(1+2*H.w.country))] <- Fti
    }
    vector.country <- vector
    
    # - extract country level
    for(i in 1:dim(smot$mts)[3]){ 
      c <- vector %*% smot$mts[,,i]
      df.mts[i, "Country"] <- c
    }
    
    
    ## for region (sum of country and region)
    for(region in regions){
      
      # - relevant vector to extract the region level
      vector <- vector.country
      vector[region] <- 1
      
      # - extract region level
      for(i in 1:dim(smot$mts)[3]){ 
        r <- vector %*% smot$mts[,,i]
        df.mts[i, region] <- r
      }
    }
    
    
    ## for farms (sum of country, region and farm (use Ftfull))
    for(i in 1:dim(smot$mts)[3]){ 
      f <- res$fullFt.list[[i]] %*% smot$mts[,,i]
      df.mts[i, rownames(f)] <- t(f)
    }
    results.list[["mts"]] <- df.mts
    
    
    # Get the smoothed variance (Cts) 
    ## create empty data frame for Cts results
    df.Cts <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + 6)) #5 regions and 1 country
    df.Cts[,1] <- D$date
    colnames(df.Cts) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
    
    ## for country
    # - relevant vector to extract the country level
    vector <- rep(0, dim(smot$Cts)[[1]])
    names(vector) <- rownames(smot$Cts[,,1])
    
    if(is.null(H.w.country)){
      vector[1] <- 1
    }else{
      Fti <- 1
      for(n in 1:H.w.country){
        Fti <- c(Fti, c(1,0))
      }
      vector[c(1:(1+2*H.w.country))] <- Fti
    }
    vector.country <- vector
    
    # - extract country level
    for(i in 1:dim(smot$Cts)[3]){ 
      c <- vector %*% smot$Cts[,,i] %*% vector
      df.Cts[i, "Country"] <- c
    }
    
    
    ## for region (considering country and region)
    for(region in regions){ 
      
      # - relevant vector to extract the region level
      vector <- vector.country
      vector[region] <- 1
      
      # - extract region level
      for(i in 1:dim(smot$Cts)[3]){ 
        r <- vector %*% smot$Cts[,,i] %*% vector 
        df.Cts[i, region] <- r
      }
    }
    
    
    ## for farms (considering country, region and farm (use Ftfull))
    for(i in 1:dim(smot$Cts)[3]){ 
      f <- res$fullFt.list[[i]] %*% smot$Cts[,,i] %*% t(res$fullFt.list[[i]]) 
      for(farm in colnames(f)){ 
        df.Cts[i, farm] <- f[farm, farm] 
        
      }
    }
    results.list[["Cts"]] <- df.Cts
  }
  return(results.list)
}



# Function to create Learning and Test sets for DLMs ----
# - df is the initial dataset
# - relevant.names are the names of the variables to be co-modeled
# - N is the time step that divides the data into the first ~3/4 of dates for Learning.set and the rest for Test.set
learning.test.sets <- function(df, relevant.names, N){
  
  # - see which nseqs != 0 are cut and move them to the closest side
  set_nseqs <- data.frame("nseq"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "start"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "end"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "set"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])))
  loop=0
  for(i in unique(df$nseq)[unique(df$nseq)!=0]){
    
    dates <- unique(df[order(df[,"date"]), "date"])
    df_nseq <- subset(df, df$nseq == i)
    min_nseq <- df_nseq[1,"date"]
    max_nseq <- df_nseq[dim(df_nseq)[1],"date"]
    
    range <- which(dates==min_nseq):which(dates==max_nseq)
    a <- intersect(range,N+1)
    
    loop=loop+1
    set_nseqs[loop, "nseq"] <- i
    set_nseqs[loop, "start"] <- which(dates==min_nseq)
    set_nseqs[loop, "end"] <- which(dates==max_nseq)
    
    if(isTRUE(length(a)>0)){
      dist_min <- N - which(dates==min_nseq)
      dist_max <- which(dates==max_nseq) - N
      
      if(dist_min < dist_max){
        set_nseqs[loop, "set"] <- "Test"
      }
      if(dist_max < dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      if(dist_max == dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      
    }else{
      
      if(isTRUE(which(dates==max_nseq) <= N)){
        set_nseqs[loop, "set"] <- "Learning"
      }else{
        set_nseqs[loop, "set"] <- "Test"
      }
    }
  }
  sum(set_nseqs$set=="Learning")/dim(set_nseqs)[1]*100   
  sum(set_nseqs$set=="Test")/dim(set_nseqs)[1]*100      
  
  learning.nseqs <- set_nseqs$nseq[set_nseqs$set== "Learning"]
  test.nseqs <- set_nseqs$nseq[set_nseqs$set== "Test"]
  
  # - for nseqs=0, use the first 3/4 of time (dates) for Learning.set
  df_0 <- subset(df, df$nseq==0)
  learning.dates.0 <- unique(df_0[order(df_0[,"date"]), "date"])[1:N]
  Learning.set.0 <- subset(df_0, df_0$date %in% learning.dates.0)
  Test.set.0 <- subset(df_0, df_0$date %in% learning.dates.0==FALSE)
  
  # - create final Learning and Test sets
  Learning.set.cut <- subset(df, df$nseq  %in% learning.nseqs)
  Learning.set <- rbind(Learning.set.0, Learning.set.cut) 
  Test.set.cut <- subset(df, df$nseq  %in% test.nseqs) 
  Test.set <- rbind(Test.set.0, Test.set.cut)  
  # - order a by site and date
  Learning.set <- Learning.set %>% 
    arrange(site, date)
  Test.set <- Test.set %>% 
    arrange(site, date)
  
  
  #remove sites in Learning.set with log0.mortality.rel.20/log0.mortality.rel.10 all equal to NA, 
  #farms with only 1 nseq (!=0) in Learning.set and 
  #all nseqs in the Learning.set with less than 6 observations != NA
  sites.exclude.all.na <- c()
  sites.exclude.1.nseq <- c()
  nseq.to.exclue.6.obs <- c()
  var <- relevant.names[length(relevant.names)]
  
  for(farm in unique(Learning.set$site)){ 
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    n.nseqs <- unique(D.farm$nseq[D.farm$nseq!=0])
    
    if(sum(is.na(D.farm[,var]))==dim(D.farm)[1]){
      a <- unique(D.farm$site)
      sites.exclude.all.na <- c(a, sites.exclude.all.na)
      
    }else{
      
      for(i.nseq in n.nseqs){ 
        nseq.farm <- subset(D.farm, D.farm$nseq == i.nseq)
        n.not.na.nseq <- dim(nseq.farm)[1] - sum(is.na(nseq.farm[,var]))
        
        if(n.not.na.nseq<6){
          nseq.to.exclue.6.obs <- c(i.nseq, nseq.to.exclue.6.obs)
        }
      }
      f.nseqs <- setdiff(n.nseqs, nseq.to.exclue.6.obs)
      if(length(f.nseqs)<2){ 
        b <- unique(D.farm$site)
        sites.exclude.1.nseq <- c(b, sites.exclude.1.nseq)
      }
    }
  }
  # - exclude nseqs with less than 6 observations != NA on the Learning.set
  Learning.set <- subset(Learning.set, Learning.set$nseq %in% nseq.to.exclue.6.obs == FALSE) 
  
  # - exclude farms in the Learning.set with log0.mortality.rel.20/log0.mortality.rel.10 always equal to NA and 
  # - farms with only 1 nseq in Learning.set  
  all.sites.exclude <- c(sites.exclude.all.na, sites.exclude.1.nseq) 
  Learning.set <- subset(Learning.set, Learning.set$site %in% all.sites.exclude == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% all.sites.exclude == FALSE)             
  
  
  #remove the farms that all environmental data is NA
  table.is.na <- c()
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    a <- c()
    relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                        "log1.d9.phyc", "log.d9.chl",
                        "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                        "log.d9.prep", "log.d9.dino", "log.d9.diato",
                        "log.d9.nano", "log.d9.pico",
                        "d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d9.no3", "log.max.daily.range.no3",
                        "log0.mortality.rel.20")
    for (name in relevant.names[1:length(relevant.names)-1]){
      
      total <- dim(D.farm)[1]
      assign(paste0("is.na.", name, sep=""), sum(is.na(D.farm[,name])))
      a <- c(a, get(paste0("is.na.", name, sep="")))
      vector <- c(farm, total, a)
    }
    table.is.na <- rbind(table.is.na, vector)
    table.is.na <- as.data.frame(table.is.na)
    colnames(table.is.na)[1:2] <- c("site", "total")
    colnames(table.is.na)[3:length(colnames(table.is.na))] <- relevant.names[1:length(relevant.names)-1]
    rownames(table.is.na) <- NULL
  }
  exclude.sites.i <- which(table.is.na$total==table.is.na$d1.temp) #if 1 environmental variable is all missing, all variables are (for this case)
  exclude.sites <- table.is.na$site[exclude.sites.i]
  Learning.set <- subset(Learning.set, Learning.set$site %in% exclude.sites== FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% exclude.sites == FALSE)            
  
  #have the same sites in Learning and Test sets
  sites.learn.not.test <- setdiff(Learning.set$site,Test.set$site) 
  sites.test.not.learn <- setdiff(Test.set$site,Learning.set$site) 
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.learn.not.test == FALSE)
  Test.set <- subset(Test.set, Test.set$site %in% sites.test.not.learn == FALSE)             
  
  return(list(Learning.set=Learning.set, Test.set=Test.set))
}



# Function to create Learning and Test sets for EM algorithm ----
# - df is the initial dataset
# - relevant.names are the names of the variables to be co-modeled
# - N is the time step that divides the data into the first ~3/4 of dates for Learning.set and the rest for Test.set
## The difference is that this one doesn't remove farms on the Learning set with only 1 nseq (!=0)
learning.test.sets.EM <- function(df, relevant.names, N){
  
  # - see which nseqs != 0 are cut and move them to the closest side
  set_nseqs <- data.frame("nseq"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "start"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "end"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "set"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])))
  loop=0
  for(i in unique(df$nseq)[unique(df$nseq)!=0]){
    
    dates <- unique(df[order(df[,"date"]), "date"])
    df_nseq <- subset(df, df$nseq == i)
    min_nseq <- df_nseq[1,"date"]
    max_nseq <- df_nseq[dim(df_nseq)[1],"date"]
    
    range <- which(dates==min_nseq):which(dates==max_nseq)
    a <- intersect(range,N+1)
    
    loop=loop+1
    set_nseqs[loop, "nseq"] <- i
    set_nseqs[loop, "start"] <- which(dates==min_nseq)
    set_nseqs[loop, "end"] <- which(dates==max_nseq)
    
    if(isTRUE(length(a)>0)){ 
      dist_min <- N - which(dates==min_nseq)
      dist_max <- which(dates==max_nseq) - N
      
      if(dist_min < dist_max){
        set_nseqs[loop, "set"] <- "Test"
      }
      if(dist_max < dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      if(dist_max == dist_min){ 
        set_nseqs[loop, "set"] <- "Learning"
      }
      
    }else{
      
      if(isTRUE(which(dates==max_nseq) <= N)){
        set_nseqs[loop, "set"] <- "Learning"
      }else{
        set_nseqs[loop, "set"] <- "Test"
      }
    }
  }
  sum(set_nseqs$set=="Learning")/dim(set_nseqs)[1]*100   
  sum(set_nseqs$set=="Test")/dim(set_nseqs)[1]*100      
  
  learning.nseqs <- set_nseqs$nseq[set_nseqs$set== "Learning"]
  test.nseqs <- set_nseqs$nseq[set_nseqs$set== "Test"]
  
  # - for nseqs=0, use the first 3/4 of time (dates) for Learning.set
  df_0 <- subset(df, df$nseq==0)
  learning.dates.0 <- unique(df_0[order(df_0[,"date"]), "date"])[1:N]
  Learning.set.0 <- subset(df_0, df_0$date %in% learning.dates.0)
  Test.set.0 <- subset(df_0, df_0$date %in% learning.dates.0==FALSE)
  
  # - create final Learning and Test sets
  Learning.set.cut <- subset(df, df$nseq  %in% learning.nseqs)
  Learning.set <- rbind(Learning.set.0, Learning.set.cut) 
  Test.set.cut <- subset(df, df$nseq  %in% test.nseqs) 
  Test.set <- rbind(Test.set.0, Test.set.cut)             
  # - order a by site and date
  Learning.set <- Learning.set %>% 
    arrange(site, date)
  Test.set <- Test.set %>% 
    arrange(site, date)
  
  
  #remove sites in Learning.set with mortality.rel.20/log0.mortality.rel.10 all equal to NA, 
  #farms with only 1 nseq (!=0) in Learning.set and 
  #all nseqs in the Learning.set with less than 6 observations != NA
  sites.exclude.all.na <- c()
  sites.exclude.1.nseq <- c()
  nseq.to.exclue.6.obs <- c()
  var <- relevant.names[length(relevant.names)]
  
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    n.nseqs <- unique(D.farm$nseq[D.farm$nseq!=0])
    
    if(sum(is.na(D.farm[,var]))==dim(D.farm)[1]){
      a <- unique(D.farm$site)
      sites.exclude.all.na <- c(a, sites.exclude.all.na)
      
    }else{
      
      for(i.nseq in n.nseqs){ 
        nseq.farm <- subset(D.farm, D.farm$nseq == i.nseq)
        n.not.na.nseq <- dim(nseq.farm)[1] - sum(is.na(nseq.farm[,var]))
        
        if(n.not.na.nseq<6){ 
          nseq.to.exclue.6.obs <- c(i.nseq, nseq.to.exclue.6.obs)
        }
      }
      f.nseqs <- setdiff(n.nseqs, nseq.to.exclue.6.obs)
      if(length(f.nseqs)<2){ 
        b <- unique(D.farm$site)
        sites.exclude.1.nseq <- c(b, sites.exclude.1.nseq)
      }
    }
  }
  # - exclude nseqs with less than 6 observations != NA on the Learning.set
  Learning.set <- subset(Learning.set, Learning.set$nseq %in% nseq.to.exclue.6.obs == FALSE) 
  
  # - exclude farms in the Learning.set with mortality.rel.20/log0.mortality.rel.10 always equal to NA
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.exclude.all.na == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% sites.exclude.all.na == FALSE)            
  
  
  #remove the farms that all environmental data is NA
  table.is.na <- c()
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    a <- c()
    relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                        "log1.d9.phyc", "log.d9.chl",
                        "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                        "log.d9.prep", "log.d9.dino", "log.d9.diato",
                        "log.d9.nano", "log.d9.pico",
                        "d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d9.no3", "log.max.daily.range.no3",
                        "log0.mortality.rel.20")
    for (name in relevant.names[1:length(relevant.names)-1]){
      
      total <- dim(D.farm)[1]
      assign(paste0("is.na.", name, sep=""), sum(is.na(D.farm[,name])))
      a <- c(a, get(paste0("is.na.", name, sep="")))
      vector <- c(farm, total, a)
    }
    table.is.na <- rbind(table.is.na, vector)
    table.is.na <- as.data.frame(table.is.na)
    colnames(table.is.na)[1:2] <- c("site", "total")
    colnames(table.is.na)[3:length(colnames(table.is.na))] <- relevant.names[1:length(relevant.names)-1]
    rownames(table.is.na) <- NULL
  }
  exclude.sites.i <- which(table.is.na$total==table.is.na$d1.temp) #if 1 environmental variable is all missing, all variables are (in this case)
  exclude.sites <- table.is.na$site[exclude.sites.i]
  Learning.set <- subset(Learning.set, Learning.set$site %in% exclude.sites== FALSE)
  Test.set <- subset(Test.set, Test.set$site %in% exclude.sites == FALSE)           
  
  #have the same sites in Learning and Test sets
  sites.learn.not.test <- setdiff(Learning.set$site,Test.set$site) 
  sites.test.not.learn <- setdiff(Test.set$site,Learning.set$site) 
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.learn.not.test == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% sites.test.not.learn == FALSE)            
  
  return(list(Learning.set=Learning.set, Test.set=Test.set))
}



## Function to get the final Learning and Test sets ----
# - df is the initial dataset
# - relevant.names are the names of the variables to be co-modeled
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - N is the time step that divides the data into the first ~3/4 of dates for Learning.set and the rest for Test.set
## Use only the sites that will be there after dividing the Learning set twice
## (so that the sites in the train.set used for the EM algorithm will be the same as in the Learning.set)
get.learning.test.sets <- function(df, relevant.names, hierarchical, N){
  
  df_0 <- df[-c(which(df$nseq==0)),]
  
  sets <- learning.test.sets(df_0, relevant.names, N)
  Learning.set_0 <- sets[["Learning.set"]]
  Test.set_0 <- sets[["Test.set"]]
  
  N.EM <- round(3*(length(unique(Learning.set_0[order(Learning.set_0[,"date"]), "date"]))/4))
  
  sets.EM <- learning.test.sets.EM(df=Learning.set_0, relevant.names, N=N.EM)
  train.set.EM_0 <- sets.EM[["Learning.set"]]
  test.set.EM_0 <- sets.EM[["Test.set"]]
  
  sites.to.remove <- setdiff(Learning.set_0$site, train.set.EM_0$site)
  sites.to.keep <- unique(train.set.EM_0$site)
  
  if(hierarchical==FALSE){
    Learning.set.final <- subset(Learning.set_0, Learning.set_0$site %in% sites.to.remove == FALSE)
    Test.set.final <- subset(Test.set_0, Test.set_0$site %in% sites.to.remove == FALSE)
  }
  if(hierarchical==TRUE){
    sets <- learning.test.sets(df, relevant.names, N)
    Learning.set <- sets[["Learning.set"]]
    Test.set <- sets[["Test.set"]]
    
    Learning.set.final <- subset(Learning.set, Learning.set$site %in% sites.to.keep == TRUE)
    Test.set.final <- subset(Test.set, Test.set$site %in% sites.to.keep == TRUE)
  }
  return(list(Learning.set=Learning.set.final, Test.set=Test.set.final))
}



##### Estimate discount factors by standard optimization  #####

# Function to be minimized (lowest RMSE) ----
# - parms is a vector with the parameters (the log transform of the observation variances and the logit transform of the discount factors)
# - vars is the number of observation variances (in this case 3 - country, region and farm)
# - D is the data set, containing the time series we want to model
# - D.months contains the information about the start and duration of the production cycles
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector with the name of the regions
getFit = function(parms, vars=3, D, D.months, mu0, C0, Spline.list, regions) {

  # Transform parameters back
  # - Variances 
  countryVar = exp(parms[1])
  regionVar = exp(parms[2])
  farmVar = exp(parms[3])

  # - Discount factors  
  deltas = c()
  for (i in (vars+1):length(parms)) {
    delta = 1.0/(1 + exp(-parms[i]))
    deltas = c(deltas, delta)
  }
  
  # Call the filtering function
  resDis = runDiscountDLM(D, D.months, mu0, C0, countryVar, regionVar, farmVar, deltas, Spline.list, regions)

  # Calculate the RMSE
  results <- extract.res(resDis, smot=NULL, D, H.w.country=NULL, regions)
  et.results <- results[["et"]]
  et.all <- na.omit(unlist(et.results[,c(2:ncol(et.results))]))
  RMSE <- round(sqrt(mean(na.omit(et.all)^2)),4)
  
  # The next line can be enabled if you wish to follow the values 
  # during optimization (but it delays the optimization)
  print(paste("Parms:", countryVar, regionVar, farmVar, toString(deltas), RMSE))
  
  return(RMSE)
}



# Function to transform the V's and deltas to log and logit, respectively ----
# - priorV is a vector of prior observation variances
# - priorDeltas is a vector of discount factors
createPriorTransform = function(priorV, priorDeltas) {
  transPrior = c()
  vars = length(priorV)
  N = vars + length(priorDeltas)
  transPrior[1:vars] = log(priorV[1:vars])
  transPrior[(vars+1):N] = log(priorDeltas/(1 - priorDeltas))
  return(transPrior)
}



# Function to transform the estimated values back to normal ----
# - optimRes is a vector with all parameters (the log transform of the observation variances and the logit transform of the discount factors)
# - vars is the number of observation variances (in this case 3 - country, region and farm)
transformResults = function(optimRes, vars) {
  pars = c()
  N = length(optimRes$par)
  pars[1:vars] = exp(optimRes$par[1:vars])
  pars[(vars+1):N] = 1.0/(1 + exp(-optimRes$par[(vars+1):N]))
  return(pars)
}



# Function to estimate the discount factors and observation variances ----
# - priorV is a vector of prior observation variances
# - priorDeltas is a vector of discount factors
# - D is the data set, containing the time series we want to model
# - D.months contains the information about the start and duration of the production cycles
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector with the name of the regions
estimateDiscountModel = function(priorVs, priorDeltas, D, D.months, mu0, C0, Spline.list, regions) {
  # How many variances
  initialParms = createPriorTransform(priorVs, priorDeltas)
  # Estimate values
  estim = optim(initialParms, getFit, gr = NULL, vars=3, D, D.months, mu0, C0, Spline.list, regions)
  # The "optim object" is returned. The parameters can be extracted by transformResults(estim, vars)
  return(estim)
}
