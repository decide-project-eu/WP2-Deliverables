#########################################################################################
# Function to calculate the cumulative sum of forecast errors for univariate models.
# The forecast errors can, optionally, be standardized.
#########################################################################################
# x: A vector of values, for which the CUSUM will be calculated
#########################################################################################
runCusumUnivariate = function(x) {
  cusum = c()
  cusum[1] = 0
  for (i in 2:length(x)) {
    cusum[i] = sum(x[2:i])
  }
  return(cusum)
}



##############################################################################
# runV.mask
# ---------
# Inputs:
#   cusums: result object returned by runCusumUnivariate(), containing cumulative sums of forecast errors.
#   dist  : V‑mask lead distance (horizontal offset of the vertex in time).
#   angle : slope of each arm of the V‑mask (in CUSUM units per time step).
#   reset : logical; if TRUE, restart CUSUM at 0 after an alarm.
##############################################################################

runV.mask <- function(cusums, dist, angle, reset = TRUE) {
  ## Pre‑allocate container -------------------------------------------------
  n    <- length(cusums)
  mask <- matrix(FALSE, nrow = n, ncol = 6)  # will store metadata
  colnames(mask) <- c("t", "CUSUM", "alarmLow", "alarmHigh", "triggerIdx", "CUSUMreset")
  
  mask[1, ] <- c(1, cusums[1], FALSE, FALSE, NA, cusums[1])
  lastAlarm <- 0            # index of most recent alarm (for reset logic)
  alarms    <- integer(0)   # vector to hold alarm time‑steps
  
  for (i in 2:n) {
    ## Populate bookkeeping columns ----------------------------------------
    mask[i, 1] <- i
    mask[i, 2] <- cusums[i]
    
    alarmPrev <- mask[i - 1, 3] || mask[i - 1, 4]  # alarm at t‑1?
    
    # Update reset CUSUM column (6) ----------------------------------------
    delta <- cusums[i] - cusums[i - 1]
    mask[i, 6] <- if (alarmPrev) delta else delta + mask[i - 1, 6]
    
    currentLevel <- if (reset) mask[i, 6] else cusums[i]
    
    ## Skip if we're immediately after a reset -----------------------------
    if (!reset || !alarmPrev) {
      for (j in (i - 1):(lastAlarm + 1)) {
        lower <- currentLevel - (i - j + dist) * angle
        upper <- currentLevel + (i - j + dist) * angle
        prev  <- if (reset) mask[j, 6] else cusums[j]
        
        if (prev < lower) {
          mask[i, 3] <- TRUE     # lower‑tail alarm
          mask[i, 5] <- j
          alarms <- c(alarms, i)
          if (reset) lastAlarm <- i
        }
        if (prev > upper) {
          mask[i, 4] <- TRUE     # upper‑tail alarm
          mask[i, 5] <- j
          if (reset) lastAlarm <- i
        }
      }
    }
  }
  
  Pred <- apply(X = mask[,c(3,4)], MARGIN = 1, FUN = max)
  
  return(list('mask' = mask, 'Pred' = Pred))
}



# A function for selecting the dist parameter for the V-mask, given the angle
choose_dist <- function(angle, h = 5) round(h / angle)


##############################################################################
# addMaskToPlot
# -------------
# Overlay a V‑mask for a single alarm point on an existing CUSUM plot.
#
# Parameters
#   mask  : matrix returned by runV.mask().
#   obs  : index (integer) at which the alarm is raised and the mask pivoted.
#   dist  : V‑mask lead distance (horizontal).
#   angle : arm slope (σ‑units per time step).
#   reset : was runV.mask() executed with reset = TRUE?
##############################################################################
addMaskToPlot <- function(mask, obs, dist, angle, reset) {
  size <- dist #obs # min(20, obs - 1)    # dist              # how many points to draw backwards
  rows <- max(nrow(mask), obs + dist)      # ensure matrix is long enough
  arms <- matrix(NA_real_, nrow = rows, ncol = 4)  # store three arm levels + x
  base <- if (reset) mask[obs, 6] else mask[obs, 2]  # CUSUM level at vertex
  
  ## Build the three arms point‑by‑point (working backwards in time) --------
  for (i in (obs + dist):(obs + dist - size)) {
    arms[i, 1] <- base - (obs + dist - i) * angle   # lower arm
    if (i > obs - 1) arms[i, 2] <- base             # central arm
    arms[i, 3] <- base + (obs + dist - i) * angle   # upper arm
    arms[i, 4] <- i                                  # x‑coordinate
  }
  
  ## Draw the arms ----------------------------------------------------------
  lines(arms[, 4], arms[, 1], lty = 2, col = "red", lwd = 1)
  lines(arms[, 4], arms[, 2], lty = 2, col = "red", lwd = 1)
  lines(arms[, 4], arms[, 3], lty = 2, col = "red", lwd = 1)
}



##############################################################################
# runAndPlotV.mask
# ----------------
# Complete wrapper: runV.mask() + plot CUSUM + overlay V‑masks for every alarm.
#
# Parameters
#   cusums: numeric vector of CUSUM values (length n >= 2)
#   dist  : V‑mask lead distance
#   angle : arm slope
#   ylims : two‑element numeric vector giving y‑axis limits for the plot
#   reset : logical; if TRUE the detector restarts after every alarm
#
# Returns
#   The mask matrix produced by runV.mask()
##############################################################################
runAndPlotV.mask <- function(cusums, dist, angle, reset = TRUE) {
  ## Run detector -----------------------------------------------------------
  Vmask.out   <- runV.mask(cusums, dist, angle, reset)  # list(mask, alarms)
  mask  <- Vmask.out$mask
  Pred <- Vmask.out$Pred
  
  ## Choose which CUSUM series to visualise --------------------------------
  y <- if (reset) mask[, 6] else mask[, 2]
  
  ## Baseline plot ----------------------------------------------------------
  plot(y, type = "b", col = "blue",
       xlab = "Time", ylab = "Cumulative sum of standardised errors",
       ylim = range(y))
  abline(h = 0)
  
  ## Overlay V‑masks only at alarm obs -----------------------------------
  alarms <- which(Pred > 0)
  for (obs in alarms) addMaskToPlot(mask, obs, dist, angle, reset)
  
  return(alarms)
}

