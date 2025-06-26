
##############################################################################
#
# A function for applying any function to a time series in a moving window
#
##############################################################################
moving.function <- function(x, n, FUN, sides=1){
  
  # Make an empty vector to store the output, which will be returned at the end
  out <- c()
  
  # Iterate over all observations in x
  for(i in 1:length(x)){
    
    # The number of sides we want to use (1 or 2) determines the window around i 
    #   we look at; this window is called j
    if(sides == 2){
      j <- (i-floor(n/2)):(i+floor(n/2))
    }else{
      j <- (i-n+1):i
    }
    
    # We can not have observations from any time before the time seriex x starts
    j <- j[which(j > 0)]
    
    # We now get the observations in the window j
    obs <- x[j]
    
    # We apply the selected function to these obserations
    res <- FUN(na.omit(obs))
    
    # We add the result to the ouput vector
    out <- c(out,res)
  }
  
  return(out)
}


##############################################################################
# Montgomery’s four Western-Electric rules with same-side checks
##############################################################################
alarmsMontgomery <- function(k, cl, SD, Ylim = NA, Main = NA, plot.it = FALSE) {
  start.time <- Sys.time()
  
  ## --- House-keeping ------------------------------------------------------
  if (is.na(Ylim[1])) {                     # sensible default y-limits
    MAX.ABS.RANGE <- max(abs(range(k)))
    Ylim <- c(-MAX.ABS.RANGE, MAX.ABS.RANGE)
  }
  if (length(SD) == 1) {                    # recycle a scalar SD
    SD <- rep(SD, length(k))
  }
  
  ## -----------------------------------------------------------------------
  ##  Rule 1 – |k – cl| > 3 SD
  ## -----------------------------------------------------------------------
  message("Applying rule 1 …")
  rule1.i      <- which(k > cl + 3*SD | k < cl - 3*SD)
  alarms_rule1 <- integer(length(k)); alarms_rule1[rule1.i] <- 1
  
  ## -----------------------------------------------------------------------
  ##  Rule 2 – Two of three successive points beyond 2 SD **on the same side**
  ## -----------------------------------------------------------------------
  message("Applying rule 2 …")
  x_above2     <- as.integer(k > cl + 2*SD)               # above CL
  x_below2     <- as.integer(k < cl - 2*SD)               # below CL
  out_above2   <- moving.function(x_above2, n = 3, FUN = sum, sides = 1)
  out_below2   <- moving.function(x_below2, n = 3, FUN = sum, sides = 1)
  rule2.i      <- which(out_above2 >= 2 | out_below2 >= 2)
  alarms_rule2 <- integer(length(k)); alarms_rule2[rule2.i] <- 1
  
  ## -----------------------------------------------------------------------
  ##  Rule 3 – Four of five successive points beyond 1 SD **on the same side**
  ## -----------------------------------------------------------------------
  message("Applying rule 3 …")
  x_above1     <- as.integer(k > cl + SD)
  x_below1     <- as.integer(k < cl - SD)
  out_above1   <- moving.function(x_above1, n = 5, FUN = sum, sides = 1)
  out_below1   <- moving.function(x_below1, n = 5, FUN = sum, sides = 1)
  rule3.i      <- which(out_above1 >= 4 | out_below1 >= 4)
  alarms_rule3 <- integer(length(k)); alarms_rule3[rule3.i] <- 1
  
  ## -----------------------------------------------------------------------
  ##  Rule 4 – Eight successive points on the same side of CL (unchanged)
  ## -----------------------------------------------------------------------
  message("Applying rule 4 …")
  x_above      <- as.integer(k > cl)
  x_below      <- as.integer(k < cl)
  out_above    <- moving.function(x_above, n = 8, FUN = sum, sides = 1)
  out_below    <- moving.function(x_below, n = 8, FUN = sum, sides = 1)
  rule4.i      <- which(out_above >= 8 | out_below >= 8)
  alarms_rule4 <- integer(length(k)); alarms_rule4[rule4.i] <- 1
  
  ## -----------------------------------------------------------------------
  ##  Combine results & (optionally) plot – everything below is as before
  ## -----------------------------------------------------------------------
  D <- data.frame(alarms_rule1, alarms_rule2, alarms_rule3, alarms_rule4)
  D$AnyRule <- apply(D, 1, max)
  
  ## ---------- P L O T T I N G  (unchanged) -------------------------------
  if(plot.it){
    layout(matrix(c(1,1,1,1,1,1,2,3), nrow = 8, byrow = FALSE))
    par(mar = c(4, 4, 1, 1))
    plot(k, type = "b", ylim = Ylim, main = Main)
    
    abline(h = cl,          col = "black",  lwd = 2)              # CL
    abline(h = cl + SD,     lty = 2, col = "blue")                # ±1 σ
    abline(h = cl - SD,     lty = 2, col = "blue")
    abline(h = cl + 2*SD,   lty = 1, col = "purple")              # ±2 σ
    abline(h = cl - 2*SD,   lty = 1, col = "purple")
    abline(h = cl + 3*SD,   lty = 2, col = "red")                 # ±3 σ
    abline(h = cl - 3*SD,   lty = 2, col = "red")
    
    rug(rule4.i, col = "green",  lwd = 3)
    rug(rule3.i, col = "blue",   lwd = 3)
    rug(rule2.i, col = "purple", lwd = 3)
    rug(rule1.i, col = "red",    lwd = 3)
    
    par(mar = c(1, 4, 1, 1))
    plot(1, type = "n", xaxt = "n", yaxt = "n",
         xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10))
    legend("center",
           legend = c("cl", "1-σ", "2-σ", "3-σ"),
           lty = c(1, 2, 1, 2), col = c("black", "blue", "purple", "red"),
           lwd = 2, horiz = TRUE, bty = "n")
    
    plot(1, type = "n", xaxt = "n", yaxt = "n",
         xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10))
    legend("center",
           legend = c("Rule 1", "Rule 2", "Rule 3", "Rule 4"),
           lty = 1, col = c("red", "purple", "blue", "green"),
           lwd = 3, horiz = TRUE, bty = "n")
    
    layout(matrix(1, 1, 1))  # reset
    
  }
  
  diff.time <- round(difftime(Sys.time(), start.time, units = "mins"), 1)
  message("Done! It took ", diff.time, " minutes.")
  
  return(D)
}


##############################################################################
#
# A function for calculating performance metrics (sensitivity, sepcificity, major mean accuracy)
# by comparing a vector of alarms with a vector of observed events 
#
# observations: the vector of observed events (binary: 1 = positive, 0 = negative)
# alarms: the vector of alarms and non-alarms (binary: 1 = alarm, 0 = non-alarm)
#
##############################################################################
getPerformance <- function(observations, alarms){
  
  # We put the two input vectors together as a data frame to make it easier to work with
  obs.n.pred <- as.data.frame(cbind('Obs'=observations, 'Pred'=alarms))
  
  # We first distinguish between the positive alarms (where the alarm is raised) 
  #  and the negative alarms (where the alarm is not raised)
  Positive.alarms <- subset(obs.n.pred, obs.n.pred$Pred == 1)
  Negative.alarms <- subset(obs.n.pred, obs.n.pred$Pred == 0)
  
  # Get the number of positve and netagive observations
  P <- nrow(Positive.alarms)
  N <- nrow(Negative.alarms)
  
  # Calculate true positive alarms (TP)
  # - this is the subset of the positive alarms where the observations are positive
  TP <- length(which(Positive.alarms$Obs == 1))
  
  # Calculate false positive alarms (FP)
  # - this is the subset of the positive alarms where the observations are negative
  FP <- length(which(Positive.alarms$Obs == 0))
  
  # Calculate true negative alarms (TN)
  # - this is the subset of the negative alarms where the observations are negative
  TN <- length(which(Negative.alarms$Obs == 0))
  
  # Calculate false negative alarms (FN)
  # - this is the subset of the negative alarms where the observations are positive
  FN <- length(which(Negative.alarms$Obs == 1))
  
  # Calculate sensitivity
  Sensitivity <- round( TP/(TP+FN) , 2)
  
  # Calculate the specificity
  Specificity <- round( TN/(TN+FP) , 2)
  
  # Calculate the major mean accuracy
  MMA <- mean(c(Sensitivity, Specificity))
  
  # Calculate the 95 % confidence intervals
  
  # standardError <- sqrt( (p * (1-p))/N  )
  
  # - for sensitivity
  standardError_Sensitivity <- sqrt( (Sensitivity * (1 - Sensitivity))/P )
  
  # - for specificity
  standardError_Specificity <- sqrt( (Specificity * (1 - Specificity))/N )
  
  # - for MMA
  standardError_MMA <- sqrt( (MMA - (1 - MMA))/(P+N) )
  
  # Calculate the 95 % CI for MMA
  lower <- round(MMA - 1.96*standardError_MMA, 3)
  upper <- round(MMA + 1.96*standardError_MMA, 3)
  CI <- paste0(lower,'-',upper)
  
  # Put it all together in an output and return it
  out <- as.data.frame(cbind('TP'=TP, 'FP'=FP, 'TN'=TN, 'FN'=FN, 'Sensitivity'=Sensitivity, 'Specificity'=Specificity, 'MMA'=MMA, 'MMA 95 % CI'=CI))
  
  return(out)
}
  

###################################################################################################################################################


#### Functions the students should make themselves ####

