
# Balance the training data 
balance_classes <- function(df, class_col) {
  # Check if the column exists
  if (!class_col %in% names(df)) {
    stop("Column not found in data frame.")
  }
  
  # Get class counts
  class_counts <- table(df[[class_col]])
  
  # Determine the target count (maximum class count)
  max_count <- max(class_counts)
  
  # List to store balanced data
  balanced_list <- list()
  
  # Perform oversampling
  for (class_name in names(class_counts)) {
    subset_df <- df[df[[class_col]] == class_name, ]
    n <- nrow(subset_df)
    
    # If the class already has the max count, keep as is
    if (n == max_count) {
      balanced_list[[class_name]] <- subset_df
    } else {
      # Sample with replacement
      additional_samples <- subset_df[sample(1:n, max_count - n, replace = TRUE), ]
      balanced_list[[class_name]] <- rbind(subset_df, additional_samples)
    }
  }
  
  # Combine all balanced class subsets
  balanced_df <- do.call(rbind, balanced_list)
  
  # Shuffle the rows
  balanced_df <- balanced_df[sample(1:nrow(balanced_df)), ]
  rownames(balanced_df) <- NULL
  
  return(balanced_df)
}

class_accuracies <- function(obs.n.pred) {
  # Check input
  if (!all(c("Obs", "Pred") %in% names(obs.n.pred))) {
    stop("Data frame must contain 'Obs' and 'Pred' columns.")
  }
  
  # Get unique classes
  classes <- sort(unique(obs.n.pred$Obs))
  
  # Store results
  accuracies <- c()

  # Compute per-class accuracy
  for (cls in classes) {
    # Subset to rows where the true class is cls
    rows <- obs.n.pred$Obs == cls
    n_total <- sum(rows)
    n_correct <- sum(obs.n.pred$Pred[rows] == cls)
    
    # Avoid division by zero
    if (n_total > 0) {
      accuracies <- c(accuracies, n_correct/n_total)
    } else {
      accuracies <- c(accuracies, NA)
    }
  }
  
  # Calculate mean accuracy (ignoring NA)
  mean_accuracy <- mean(na.omit(accuracies))
  
  # Combine into a named vector
  result <- c(accuracies, Mean = mean_accuracy)
  result <- as.data.frame(t(as.data.frame(result)))
  colnames(result) <- c(classes, 'Mean')
  rownames(result) <- NULL
  return(result)
}


#### A function for seeing the seperation of the classes achieved by the trained model
see.seperation <- function(obs.n.pred_orig, Main=NA, ylim=c(0,1)){
  
  # Plot a histogram showing the distribution of predictions given the true state
  # plot bar chart for age with Q1-Q3 intervals around the median
  summaries <- as.data.frame(as.data.frame(aggregate(obs.n.pred_orig$Pred, by=list(obs.n.pred_orig$Obs), FUN=summary)))
  labels <- summaries$Group.1
  values <- as.data.frame(summaries$x)
  Medians <- values$Median
  lower <- values$`1st Qu.`
  upper <- values$`3rd Qu.`
  
  Medians <- rbind(Medians)
  colnames(Medians) <- labels
  BP <- barplot(Medians, xaxt = "n", beside = TRUE,  ylim=ylim, col=c( 'grey80', 'grey30'), space=c(0.2,0), las=2, main=Main)
  # - ad rotated labels
  arrows(x0=BP, y0=lower, x1 = BP, y1 = upper, code = 3, angle = 90, length = 0.10)
  text(x = 1:length(Medians)-0.5,
       y = -0.05,
       labels = colnames(Medians),
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 0,
       adj = 1)
  
  # Get the range of the medians as a proxy for how well the groups are seperated
  
  return(diff(range(Medians)))
  
}

# Transform the obs.n.pred data into something that can be used by ggplot
# Gathering columns a to j
confusion.matrix <- function(obs.n.pred, Title=NA){
  
  library("ggplot2")
  library("dplyr")
  
  PT <- cbind()
  for(obs in unique(obs.n.pred$Obs)){
    for(pred in unique(obs.n.pred$Pred)){
      obs.pred.set <- subset(obs.n.pred, obs.n.pred$Obs == obs & obs.n.pred$Pred == pred)
      freq <- nrow(obs.pred.set)/length(which(obs.n.pred$Obs == obs))
      PT <- rbind(PT, cbind('Observation'=obs, 'Prediction'=pred, 'Freq'=freq))
    }
  }
  PT <- as.data.frame(PT)
  PT$Freq <- round(as.numeric(as.character(PT$Freq)),2)
  
  # This next step is to ensure that the plot has the accuracies in the correct diagonal
  # - (from upper left corner to lover right corner)
  
  PT <- PT %>%
    mutate(Observation = factor(Observation, levels = sort(unique(PT$Observation)) ), # alphabetical order 
           Prediction = factor(Prediction, levels = rev(sort(unique(PT$Prediction))) )) # force reverse alphabetical order
  # - this way, the diagonal represents the accuracy of the prediction of each class,
  #   and the off-diagonal values represent what proportion of each observation of each class was predicted as each of the other classes
  
  
  # Now make the fancy plot
  myplot <- ggplot(PT, aes(x=Observation, y=Prediction, fill=Freq)) +
    geom_tile() + theme_bw() + coord_equal() +
    ggtitle(Title) + 
    scale_fill_distiller(palette="Blues", direction=1) +
    guides(fill=F) + # removing legend for `fill`
    geom_text(aes(label=Freq), color="black") # printing values
  
  print(myplot)
  
  
}


# Make a nice plot to show the effect of the different strategies
cld.plot <- function(model_means, ylim=c(0,1)){
  
  # add letters to each mean
  CLD <- cld(object = model_means,
             adjust = "sidak",
             Letters = letters,
             alpha = 0.05)
  
  # Plot it
  old.par <- par(mai=c(1,1,1.25,1), no.readonly=TRUE)
  plot(CLD, horizontal = FALSE, colors = c('red', 'blue'), ylab = NULL, cex=2, ylim=ylim)
  
  # return(CLD)
}





# A function to give you the number of unique elements in a vector
get.N.unique <- function(x){
  n <- length(unique(x))
  return(n)
}

# A function for assigning observations in a data set to one of N folds for cross-validation
assign_Folds <- function(D, stratify.by, N_folds, identifyer){
  
  N.unique <- length(unique(D[,stratify.by]))
  N.per.fold <- floor(N.unique/N_folds)
  Fold <- rep(NA, nrow(D))
  D <- cbind(Fold, D)
  for(i in 1:N_folds){
    agg <- aggregate(D$Fold, by=list(D[,identifyer]), FUN=mean)
    na.i <- which(is.na(agg[,2]))
    na.obs <- agg[na.i, 1]
    if(i == N_folds){
      fold.obs.i <- which(D[,stratify.by] %in% na.obs)
      D$Fold[fold.obs.i] <- i
    }else{
      fold.obs <- sample(x = na.obs, size = N.per.fold, replace = FALSE)
      fold.obs.i <- which(D[,stratify.by] %in% fold.obs)
      D$Fold[fold.obs.i] <- i
    }
  }
  
  agg <- aggregate(D[,stratify.by], by=list(D[,'Fold']), FUN=get.N.unique)
  colnames(agg) <- c('Fold', paste0('Unique_', stratify.by))
  print(agg)
  
  return(D)
  
}


# A function for applyin the roc function to a set of observations and corresponding probability predictions extracting relevant results
apply.roc <- function(obs.test, pred.test){
  
  # Get the ROC cuve of RF applied to the test set, and extract the auc
  roc.out.test <- roc(response = obs.test, predictor = pred.test)
  AUC <- round( roc.out.test$auc , 2)
  
  # Calculate 95 %  CI confidence interval
  
  # - get the number of positive observations
  N1 <- length(which(obs.test == 1))
  
  # - get the number of negative observations
  N2 <- length(which(obs.test == 0))
  
  # - do the rest
  Q1 <- AUC/(2-AUC)
  Q2 <- (2*AUC^2)/(1+AUC)
  SE.AUC <- sqrt( ( AUC*(1-AUC)+(N1-1)*(Q1-AUC^2)+(N2-1)*(Q2-AUC^2) ) /(N1*N2) )
  
  AUC.min <- round(AUC - 1.96*SE.AUC,2)
  AUC.max <- round(AUC + 1.96*SE.AUC,2)
  AUC.CI <- paste(AUC.min, '-', AUC.max, sep='')
  
  # Get the three vecotors "thresholds", "sensitivities", and "specificities"
  thresholds <- round( roc.out.test$thresholds , 3)
  sensitivities <- round( roc.out.test$sensitivities , 2)
  specificities <- round( roc.out.test$specificities , 2)
  
  # Calculate the MMA for all thresholds
  MMA.all <- (sensitivities+specificities)/2 
  
  # Find the maximum MMA value
  MMA.max <- max(MMA.all)
  
  #Find the threshold which results in the highest MMA
  best.threshold <- thresholds[which(MMA.all == MMA.max)] 
 
  # If we have more than one threshold which result in the best MMA, we assume that a higher specificity is better than a higher sensitivity!
  best.threshold <- best.threshold[length(best.threshold)]
  best.i <- which(thresholds %in% best.threshold)
  
  # Find the sensitivity and specificity associated with the maximum MMA
  Se.best <- sensitivities[best.i]
  Sp.best <- specificities[best.i]
  
  # Plot the ROC curve
  plot(roc.out.test )
  abline(h=Se.best, lty=2)
  abline(v=Sp.best, lty=2)
  
  # Look 
  out <- as.data.frame(cbind('AUC'=AUC, 'AUC.CI'=AUC.CI, 'MMA.max'=MMA.max, 'Se'=Se.best, 'Sp'=Sp.best, 'Threshold'=best.threshold ))
  return(out)
  
}
