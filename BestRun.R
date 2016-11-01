## Added diff pattern solve and changed diffTableSolve to add ability to return 0 pts of prediction
require(e1071)
require(numbers)
set.seed(123)

## Params
# Min points to attempt fit
mp <- 2
# pc deviation allowed before reverting to mode
pcThresh <- 0.01
pcThreshSVM <- 0.02
# Degrees of fits to try
# degs <- 1:6
degs <- 0 # off
# n previous pts to try
# pps <- 2:3
pps <- 1:12 # off
# Plot when fitting? May not always work
plotOn <- FALSE
# SVM options
ppSVMOn <- FALSE
degSVMOn <- FALSE
diffOn = TRUE
acceptDodgy = TRUE
pattOn = TRUE
patMatchThresh = 0.95 # Internal  (1 good)
# patMatchThresh2 = 0.0001 # Extrernal (RRScore: 0 good) Not added yet
# Recurrance relation solve
RROn = TRUE
# Required RR score to use (from score PC)
RRScore = 0 # External
# Difftable pattern
diffTablePatternOn = TRUE
diffTablePatternThresh = 0.0001 # External

## Functions 

RRSolve <- function(y, n)
{
  # seq <- c(1,5,11,21,39,73,139,269,527)
  # n <- 3L
  i <- rep(1L:n, n) + rep(1L:n, each = n) - 1L
  A <- matrix(y[i], nrow = n, ncol = n)
  j <- (n + 1L):(2L * n)
  B <- matrix(y[j], nrow = n, ncol = 1L)
  # A * x = B
  
  # x <- solve(A, B)
  
  x <- tryCatch({
    x <- solve(A, B)
  }, error = function(err){
    x <- NaN
    return(x)
  }, finally = {
  })
  
  if (!is.na(x)[1])
  {
    terms <- tail(y, n)
    # Alternate code:
    # prediction <- crossprod(x, terms)[[1]]
    prediction <- sum(x  * terms) 
    # Answer 1041
    return(list(c(y, prediction), prediction))
  } else
    return(list(c(y, NaN), NaN))
}

patternSearch <- function(y, thresh=0.98)
{
  m = length(y)
  for (p in 1:floor(m/2)) 
  {
    # n full multiples of p
    nM = floor(m/p)
    # Get end of sequence
    yCheck = y[((m-nM*p)+1):m]
    pCheck = y[((m-p)+1):m]
    
    l = yCheck==pCheck
    
    propMatched = sum(l)/length(l)
    
    # print(yCheck)
    # print(pCheck)
    # print(l)
    #print(propMatched)
    if (propMatched > thresh)
    { 
      #   cat("Success \n")
      return(list(TRUE,pCheck,propMatched))
    }
  }
  #cat("Fail")
  return(list(FALSE,NaN,0))
}

# nPts = 0 won't add any predictions
diffTableSolve2 <- function(y, depth=1, nPts=1) 
{
  # Attempt to use difference table to solve sequence
  yPred <- NaN
  d = list()
  d[[1]] = y
  i = 1;
  conv <- FALSE
  convFail <- FALSE
  goodness <- "??"
  m=length(y)
  while ((!conv & !convFail))
  {
    i = i+1
    # cat(i)
    # thisDiff <-diff(d[[i-1]])
    
    # Calculate this diff rolling window
    dm = m-(i-1)
    thisDiff <- matrix(nrow = 1, ncol = dm)
    for (de in 1:dm)
    {
      fp = de
      np = de+depth
      thisDiff[de] <- d[[i-1]][np]-d[[i-1]][fp]
      
    }
    
    thisDiff<- thisDiff[!is.na(thisDiff)]
    d[i] <-list(thisDiff)
    
    # Check for convergence
    conv = length(unique(d[[i]]))==1 & length(thisDiff)>1
    if (conv) {
      if (length(d[[i]])<=(2*depth))
      {
        goodness = "Dodgy"
      } else {
        goodness = "OK"
      }
    } else { # Check for failure
      if (i == length(y) | length(d[[i]])==1 | length(d[[i]])<(depth))
      {
        convFail = TRUE
        goodness = "Failed"
      }
    }
    
    
    ## If at end, fail
    #if (i == length(y) | length(d[[i]])==1 | length(d[[i]])<(depth/1.1))
    #{
    #  convFail = TRUE
    #  goodness = "Failed"
    #} else # Check for convergence
    #{
    #  conv = length(unique(d[[i]]))==1
    #  if (length(d[[i]])<=(2*depth))
    #  {
    #    goodness = "Dodgy"
    #  } else {
    #    goodness = "OK"
    #  }
    #}
  }
  
  d = (list(conv, yPred, d, depth, goodness)) 
  if (nPts>0)
    # Add preditiction if requested
  {
    d = diffTablePredict2(d, depth, nPts)
  }
  return(d)
}

diffTablePredict2 <- function(d,depth,nPts)
{
  conv = d[[1]]
  yPred = NaN
  i = length(d[[3]])
  m = length(d[[3]][[1]])
  d2 = d[[3]]
  for (dd in 1:nPts)
  {
    d2[[i]][length(d2[[i]])+1]=d2[[i]][1]
    if (conv) # Converged
    {
      # Work back
      for (p in i:2)
      {
        d2[[p-1]][length(d2[[p-1]])+1] = d2[[p]][length(d2[[p]])] + d2[[p-1]][length(d2[[p-1]])-(depth-1)]
      }
      # Add last
      yPred = d2[[1]][(m+1):(m+nPts)]
    }
  }
  
  d[[3]] = d2
  d[[2]] = yPred
  return(d)
}

diffTablePattern <- function(y, mThresh = 0.97, mPatLen2Match = 4)
{
  # Rerun diffPattern solve at all depths - ignore convergence
  # Run pattern solve on all levels - record scores (threshold = 0)
  # Find best score (at lowest depth and level): Better than mThresh?
  # Redo difftable solve and pattern search on best
  # Calculate next value
  m = length(y)
  nDepths = m
  
  if (m<=mPatLen2Match)
  {
    # Too short, don't run
    result = list(FALSE, NaN, NaN, "Too short", NaN)
  } else {
    
    # Prepare matrix to save patterns
    # Won't use all rows a n patterns reduces as depth increases
    pMat = matrix(nrow = m, ncol = m)
    
    for (d in 1:nDepths)
    {
      DTS = diffTableSolve2(y, d, nPts=0)
      
      pats = DTS[[3]]
      nPats = 1:length(pats)
      for (p in 1:length(nPats))
      {
        # One is the actual sequence, could exclude to optimise
        pat = DTS[[3]][[p]]
        
        # Only score if final level is long enough (not "dodgy")
        if (length(pat) >= mPatLen2Match)
        {
          thisScore = patternSearch(pat, 0)[[3]]
          pMat[d,p] = thisScore
        }
      }
    }
    
    # Find the best value to use
    # Ignore first columns because that's just the sequence at the moment
    pMat[,1] = -1
    # Best score
    bSc =  max(pMat, na.rm=TRUE)
    # Good enough?
    if (bSc >= mThresh)
    {
      # Where?
      coords = which(pMat == bSc, arr.ind = TRUE)
      d = coords[[1,1]]
      p = coords[[1,2]]
      
      # Redo calculations
      DTS = diffTableSolve2(y, d, nPts=0)
      pat = DTS[[3]][[p]]
      PTS = patternSearch(pat, 0)
      pat2Add = PTS[[2]]
      
      # Add the pattern to its level
      lm = length(DTS[[3]][[p]])
      pm = length(pat2Add)
      # DTS[[3]][[p]][(lm+1):(lm+1+pm)] = pat2Add
      # Or just use concatonate...
      DTS[[3]][[p]] = c(DTS[[3]][[p]], pat2Add)
      
      # Now work back and get preditcion
      # Doing manually rather than using diffTablePredict2
      
      # Start from level the pattern was found on
      # Ditch everything else
      newTable = DTS[[3]][1:p]
      nv = -1
      for (l in (p-1):1)
      {
        nv = nv+1
        # On each level, just add one more value for now
        # So even if longer pattern was added, just take the previous end value (remember depth!)
        # Length of this level:
        tlm = length(newTable[[l]])
        newTable[[l]] = c(newTable[[l]], newTable[[l]][(tlm+1)-d] + newTable[[l+1]][lm+nv])
      }
      
      # Should be n(y) + 1
      newYLen = length(newTable[[1]])
      yPred = newTable[[1]][newYLen]
      result = list(TRUE, newTable, yPred, pMat, coords)
      
    } else {
      # Failed
      
      result = list(FALSE, NaN, NaN, pMat, NaN)
    }
  }
  
  return(result)
}

# Mode function
# http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
mode <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Find greatest common divisor in sequnece
# Too slow, use mGCD in numbers instead
findGCD <- function(y)
{
  sMax = max(y)
  gcd <- NaN
  g = 0
  nInt = sMax
  stop <- FALSE
  while (!stop)
  {
    cat(nInt, "\n")
    if (all((y/nInt) == round(y/nInt)))
    {
      gcd <- nInt
      stop <- TRUE
    } else {
      g <- g+1
      nInt = sMax-g
      if (nInt==0) {stop<-TRUE}
    }
  }
  return(gcd)
}

# Calculate sequence sig
seqSig <- function(y)
{
  if (y[y!=0][1]>0) {sign=1} else {sign=-1}
  
  # gcd <- findGCD(y)
  gcd <- mGCD(y)
  
  sig = sign*y/gcd
  
  return(list(sig, gcd, sign))
  
}

diffTableSolve <- function(y) 
{
  # Attempt to use difference table to solve sequence
  yPred <- NaN
  d = list()
  d[[1]] = y
  i = 1;
  conv <- FALSE
  convFail <- FALSE
  while ((!conv & !convFail))
  {
    i = i+1
    # cat(i)
    thisDiff <-diff(d[[i-1]])
    
    d[i] <-list(thisDiff)
    
    
    # If at end, fail
    if (i == length(y))
    {
      convFail = TRUE
    } else # Check for convergence
    {
      conv = length(unique(d[[i]]))==1
    }
  }
  
  if (conv) # Converged
  {
    # Work back
    for (p in i:2)
    {
      d[[p-1]][length(d[[p-1]])+1] = d[[p]][length(d[[p]])] + d[[p-1]][length(d[[p-1]])]
    }
    # Add last
    yPred = d[[1]][length(d[[1]])]
    # Don't bother scoring
    # Just assume if it convereged, it's correct
    
  }
  return(list(conv, yPred, d))  
}

# Score PC for use with fits (manages nans) -but why?
scorePC <- function(yP, yOr)
{
  # Closer to 0 is better
  score = abs(1-yP/yOr)
  
  # Note 0/0 = nan, 1/0 = inf
  if (is.infinite(score)) {score=999}
  # if (is.nan(score)) {score=0} # What was the logic here?
  if (is.nan(score)) {score=8989}
  
  return(score)
}

scorePCRR <- function(yP, yOr)
{
  # Closer to 0 is better
  score = abs(1-yP/yOr)
  
  # Note 0/0 = nan, 1/0 = inf
  if (is.infinite(score)) {score=999}
  # if (is.nan(score)) {score=0}
  
  return(score)
}

scoreMSE <- function(yP, yOr)
{
  # Closer to 0 is better
  score = 1/length(yP) * sum(((yP-yOr)^2))
  
  # Scale
  score = score/(max(yOr)-min(yOr))
  
  # Note 0/0 = nan, 1/0 = inf
  # if (is.infinite(score)) {score=999}
  # if (is.nan(score)) {score=0}
  
  return(score)
}

expandFeaturesNL <- function(xy, fStr) 
{
  # If y is in df, remove
  holdY = ('y' %in% names(xy))
  if (holdY)
  {
    y <- xy['y']
    xy$y <-NULL
  }
  
  # Check there are enough columns
  nfs = ncol(xy)
  if (nfs>1) 
  {
    # For each column...
    for (c1 in 1:(nfs))
    {
      # * and / with all other columns and update names
      for (c2 in 1:nfs)
      {
        if (!c1==c2)
        {
          # Add this column times the other columns
          cn = paste(colnames(xy)[c1], 'm', colnames(xy)[c2], sep="")
          fStr <- paste(fStr, "+", cn, sep="")
          #need to update to get correct colu,n using name
          xy <- data.frame(xy, xy[,c1]*xy[,c2])
          colnames(xy)[ncol(xy)] <- cn
          
          # Add this column divided by the other columns
          cn = paste(colnames(xy)[c1], 'd', colnames(xy)[c2], sep="")
          fStr <- paste(fStr, "+", cn, sep="")
          #need to update
          xy <- data.frame(xy, xy[,c1]/xy[,c2])
          colnames(xy)[ncol(xy)] <- cn
        }
      }
    }
  }
  
  # Check
  xy[is.na(xy)] = 0
  xy[xy==Inf] = 0
  xy[xy==-Inf] = 0
  
  # If y was removed, return
  if (holdY)
  {
    xy <- data.frame(y, xy)
  }
  
  return(list(xy, fStr))
}  

expandFeaturesLin <- function(xy)
{
  nfs = ncols(xy)-1
  for (c in 2:(nfs+1))
  {
    # Add sin, cos, log, exp etc
  }
}

# Fit, tune, and score SVM with m points
fitSVMDegs <- function(y, m, plotOn=FALSE)
{
  # Fits SVM as degree of x
  # Fits SVM (skip to save time?)
  # Roughly tunes SVM - loss function is MSE
  # Retuned with higher resolution - loss function is still MSE
  # Use best model to predict yPred
  # Scores this model with scorePC, and returns next prediction
  yOr <- sequ
  y <- sequ[1:m]
  x <- 1:m
  
  df = data.frame(X=x, Y=y)
  df2 = df
  df2[m+1,] <- c(m+1, NaN)
  
  # f <- svm(Y ~ X , df)
  # 
  # yPred <- predict(f, df2)
  # 
  # if (plotOn)
  # {
  #   gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
  #   plot(df2$X, yPred, col = "red", pch=4, ylim=gYLim)
  #   points(x,y)
  #   
  # }
  
  # perform a grid search
  fTuned <- tune(svm, Y ~ X,  data = df, kernal="polynomial",
                 ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:10)))
  
  if (plotOn) {
    plot(fTuned)
    yPred <- predict(fTuned$best.model, df2)
    gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
    plot(df2$X, yPred, col = "red", pch=4, ylim=gYLim)
    points(x,y)
  }
  
  eps = fTuned$best.model$epsilon
  lowerEps = eps-0.1;
  upperEps = eps+0.1;
  # Limit
  if (lowerEps<0) {lowerEps=0}
  if (upperEps>1) {upperEps=1}
  
  cost = fTuned$best.model$cost
  upperCost = cost+20
  lowerCost = cost-20
  # Limit
  if (lowerCost<0) {lowerCost=1}
  
  fTuned <- tune(svm, Y ~ X,  data = df, kernal="polynomial",
                 ranges = list(epsilon = seq(lowerEps,upperEps,0.05), cost = seq(lowerCost, upperCost,5)))
  
  if (plotOn) {
    plot(fTuned)
    yPred <- predict(fTuned$best.model, df2)
    gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
    plot(df2$X, yPred, col = "red", pch=4, ylim=gYLim)
    points(x,y)
  }
  
  # Get best model
  f <- fTuned$best.model
  
  # Add next point to predict
  df = rbind(df, c(m+1,NaN))
  yPred <- predict(f, df) 
  
  if (plotOn)
  {
    gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
    plot(df$X,yPred, col="red", ylim=gYLim)
    points(x,y)
    
  }
  
  
  score <- scorePC(yPred[m], yOr[m])
  
  return(list(f, score, yOr, round(yPred)))
}

fitSVMPP <- function(y, m, pp, plotOn=FALSE)
{
  yOr <- sequ
  y <- sequ[1:m]
  x <- 1:m
  
  # Ditch start of y, where not enough previous points are available
  y <- data.frame(y=yOr[(pp+1):m])
  x <- (pp+1):m
  m2 <- length(x)
  
  # Create forumula and data frame for fitting
  fStr <- "y~"
  for (ppi in 1:pp)
  {
    # For each previous value of y
    # Add to formulaString
    fStr <- paste(fStr, "+x", ppi, sep="")
    
    # Create new column (rolling window of previous values)
    nC <- yOr[(pp-(ppi-1)):(m-ppi)]
    # Append to data frame
    y <- data.frame(y, nC)
    # Apply column name
    colnames(y)[ppi+1] <- paste("x", ppi, sep ="")
  }
  
  df = y
  df2 = y
  # df = data.frame(X=x, Y=y)
  # f <- svm(y ~ x1 , df)
  # fstr2 <-y ~ x1+x2
  
  f <-svm(y ~., y, kernal="polynomial", degree=pp)
  df2[m2+1,] <- c(NaN, yOr[m:((m+1)-(pp))])
  yPred <- predict(f, newdata=df2)
  
  if (plotOn)
  {
    gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
    plot((pp+1):(m+1), yPred, col = "red", pch=4, ylim= gYLim)
    points(x,df$y)
    
  }
  
  # perform a grid search
  fTuned <- tune(svm, y ~ .,  data = df, kernal="polynomial", degree=pp,
                 ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:10)))
  
  #gamma = 2^(-1:1),
  
  
  if (plotOn) {
    gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
    plot(fTuned, ylim= gYLim)
    yPred <- predict(f, df2)
  }
  
  eps = fTuned$best.model$epsilon
  lowerEps = eps-0.1;
  upperEps = eps+0.1;
  # Limit
  if (lowerEps<0) {lowerEps=0}
  if (upperEps>1) {upperEps=1}
  
  cost = fTuned$best.model$cost
  upperCost = cost+20
  lowerCost = cost-20
  # Limit
  if (lowerCost<0) {lowerCost=1}
  
  fTuned2 <- tune(svm, y ~ .,  data = df,kernal="polynomial", degree=pp,
                  ranges = list(epsilon = seq(lowerEps,upperEps,0.02), cost = seq(lowerCost, upperCost,4)))
  
  if (plotOn) {plot(fTuned)}
  
  # Get best model
  f <- fTuned2$best.model
  
  # Add next point to predict
  df[m2+1,] <- c(NaN, yOr[m:((m+1)-(pp))])
  
  yPred <- predict(f, df) 
  
  if (plotOn)
  {
    gYLim = range(min(min(yOr),min(yPred)),max(max(yOr),max(yPred)))
    plot((pp+1):(m+1),yPred, col="red", ylim= gYLim)
    points(x,yOr[(pp+1):x[length(x)]], pch=5)
  }
  
  
  score <- scorePC(yPred[m-pp], yOr[m])
  
  return(list(f, score, yOr, round(yPred)))
}

# Fit and score model using previous y function
fitModPP <- function(y, m, pp, expFea=FALSE, plotOn=FALSE)
{
  # Fit using m points of y with polynomial degree d
  # if m<length(y) uses next known point to score fit
  yOr <- y
  # Create data frame to use with fit()
  # Ditch start of y, where not enough previous points are available
  y <- data.frame(y=yOr[(pp+1):m])
  x <- (pp+1):m
  m2 <- length(x)
  
  # Create forumula and data frame for fitting
  fStr <- "y~"
  for (ppi in 1:pp)
  {
    # For each previous value of y
    # Add to formulaString
    fStr <- paste(fStr, "+x", ppi, sep="")
    
    # Create new column (rolling window of previous values)
    nC <- yOr[(pp-(ppi-1)):(m-ppi)]
    # Append to data frame
    y <- data.frame(y, nC)
    # Apply column name
    colnames(y)[ppi+1] <- paste("x", ppi, sep ="")
  }
  
  if (expFea)
  {
    expY = expandFeaturesNL(y, fStr)
    y <- expY[[1]]
    fStr <- expY[[2]]
  }
  
  # Do fit
  f <- lm(formula(fStr), data=y)
  
  # Predict next point and score fit, if possible
  newData = as.data.frame(y[,2:(pp+1)])
  colnames(newData) <- colnames(y)[2:(pp+1)]
  newData[m2+1,] <- yOr[m:((m+1)-(pp))]
  if (expFea)
  {
    newData <- expandFeaturesNL(newData, fStr)
    newData <- newData[[1]]
  }
  
  yP <- round(predict(f, newdata=newData))
  
  # Score
  score <- scorePC(yP[m2+1], yOr[m+1])
  # score <- scoreMSE(yP, yOr[(length(yOr)-m2):length(yOr)])
  
  
  if (plotOn) # Needs updating to work
  {
    plot(1:(m+1), yOr)
    
    lines(1:(m+1), yP)
  }
  
  # Return fit, score, y (in data frame) in a list
  return(list(f, score, y, yP))
}

# Fit and score model function
fitModDeg <- function(y, m, d, expFea=FALSE, plotOn=FALSE)
{
  # Fit using m points of y with polynomial degree d
  # if m<length(y) uses next known point to score fit
  
  yOr <- y
  # Create data frame to use with fit()
  y <- data.frame(y=yOr[1:m])
  x = 1:m
  
  # Create forumula 
  fStr <- "y~"
  for (di in 1:d)
  {
    # For each degree of x
    # Add to formulaString
    fStr <- paste(fStr, "+x", di, sep="")
    # Append column of x^deg to data frame
    y <- data.frame(y, x^di)
    # And name the column
    colnames(y)[di+1] <- paste("x", di, sep ="")
  }
  
  if (expFea)
  {
    expY = expandFeaturesNL(y, fStr)
    y <- expY[[1]]
    fStr <- expY[[2]]
  }
  
  # Do fit
  f <- lm(formula(fStr), data=y)
  
  # Predict next point and score fit, if possible
  newData = as.data.frame(y[,2:(d+1)])
  # Make sure colnames are copied (if deg=1, above gives matrix hence as.data.frame - this drops col names)
  colnames(newData) <- colnames(y)[2:(d+1)]
  newData[m+1,] <- (m+1)^degs[1:d]
  if (expFea)
  {
    newData = expandFeaturesNL(newData, fStr)
    newData <- newData[[1]]
  }
  
  yP <- round(predict(f, newdata=newData))
  
  # Score
  score <- scorePC(yP[m+1], yOr[m+1])
  #score <- scoreMSE(yP, yOr)
  
  if (plotOn)
  {
    plot(1:(m+1), yOr)
    
    lines(1:(m+1), yP)
  }
  
  # Return fit, score, y (in data frame) in a list
  return(list(f, score, y, yP))
  
}

######################################################################## Run
# Import
data <- read.csv("test.csv", stringsAsFactors=FALSE)

n <- nrow(data)
# n = 100
yPred <- matrix(nrow=n, ncol=5)
# For each sequence
for (r in seq(1,n))
{
  
  # r = 55423 - interesting sequence  
  # r = 301 - false convergence sequence?
  # r = 84099 two-state sequence
  # r = 58645 abs(sine) sequence + noise
  # r = 543 solved with Diff pattern solve, or at depth 32 (but "dodgy")
  
  # Get sequence and convert to numeric
  sequ <- data[r,2]
  sequ <- as.numeric(strsplit(sequ, ",")[[1]])
  m <- length(sequ)
  s <- ""
  
  # Blank these to be safe
  useTableVal <- FALSE
  useDegsFit <- FALSE 
  usePPSFit <- FALSE
  useFit <- FALSE
  useExp <- FALSE
  useSVM <- FALSE
  useTableVal2 = FALSE
  useTableVal3 = FALSE
  usePattern = FALSE
  useDiffPattern <- FALSE
  bestScore <- NaN
  dBestScore <- NaN
  svmScore = NaN
  ppBestScore <- NaN
  useRRVal <- FALSE
  useDTP = FALSE
  pc <- 0
  stop <- FALSE
  res <- 0
  resPat <- 0
  s <- ""
  
  wwd = NaN
  # 1212: RR
  # 6565: DTP
  # 1: diffTableSolve
  # 2: diffTableSolve on seq sig
  # 30: Difftable solve with extra degree
  # 31: Difftable solve with extra degree (dodgy)
  # 40: Pattern search
  # 3: Deg fit
  # 4: Deg fit with NL features
  # 5: PP fit
  # 6: PP fit with NL features
  # 7: Deg SVM
  # 8: PP (1) SVM
  
  # Skip straight to mode?
  if (m < mp)
  {
    # Too few points, just use mode
    stop <- TRUE
  }   
  
  # First try doing table solve
  if (!stop & diffOn)
  {
    # sequ = c(9,73,241,561,1081,1849)
    res = diffTableSolve2(sequ,1)
    
    if (res[[1]]) 
    {
      useTableVal <- TRUE
      stop <- TRUE
      wwd = 1
    }
  }
  
  # Next try after take sequence sig
  if (!stop & diffOn)
  {
    # sequ = c(9,73,241,561,1081,1849)
    ss = seqSig(sequ)
    gcd = ss[[2]]
    # Only continue if not just repeating first try
    if (!is.nan(ss[[2]]) & !gcd==1)
    {
      res = diffTableSolve2(ss[[1]],1)
      if (res[[1]]) 
      {
        # Convert result back
        res[[2]] = res[[2]] * gcd * ss[[3]]
        useTableVal <- TRUE
        stop <- TRUE
        wwd = 2
      }
    }
  }
  
  # Try difftable solve at all degrees
  dep=0
  stopDiff <- FALSE
  if (!stop & diffOn)
  {
    while (!stopDiff)
    {
      dep=dep+1
      if (dep>(m-2))
      {
        stopDiff = TRUE
        useTableVal = FALSE
      } else {
        # cat(dep, "\n")
        a = diffTableSolve2(sequ, dep)
        
        if (a[[1]] & a[[5]]=="OK")
        {
          useTableVal2 = TRUE
          useTableVal3 = FALSE
          #s <- paste(s, "Diff table solved with degree", dep, "(", a[[5]], ")", sep=" ")
          stopDiff <-TRUE
          stop <- TRUE # Good found, stop futher tests
          res=a
          wwd = 30
        } else if (a[[1]] & a[[5]]=="Dodgy" & is.na(yPred[r,1]))
        { # Mark answer if one hasn't been entered (prefer lower depth), but don't stop
          useTableVal3 = TRUE
          useTableVal2 = FALSE
          #s <- paste(s, "Diff table solved with degree", dep, "(", a[[5]], ")", sep=" ")
          stopDiff <- FALSE
          # Good not yet found, keep going
          res = a
          wwd = 31
        }
      }
    }
  }
  ## UPDATE ABOVE - is res the same as before? Might be more than one digit? Update the useTableVal switch and it's behaviour below
  
  # Pattern search
  if (!stop & pattOn)
  {
    resPat = patternSearch(sequ, thresh = patMatchThresh)
    if (resPat[[1]])
    {
      # Pattern found in sequence above threshold
      usePattern = TRUE
      useTableVal3 = FALSE
      stop <- TRUE
      
      wwd = 40
    }
  }
  
  # Run diffTable Pattern
  if (!stop & diffTablePatternOn)
  {
    resDTP = diffTablePattern(sequ[1:(m-1)])
    sc = scorePCRR(sequ[m], resDTP[[3]][1])
    if (!is.na(sc) & sc<=diffTablePatternThresh)
    {
      # Now make sure pattern isn't voided by adding 1 more point...
      
      resDTP = diffTablePattern(sequ)
      sc = scorePCRR(sequ[m], resDTP[[3]][1])
      if (!is.na(sc) & sc<=diffTablePatternThresh)
      {
        useDTP = TRUE
        wwd = 16
        stop = TRUE
      } else {
        useDTP = FALSE
      }
    } else {
      useDTP = FALSE
    }
  }
  # If a "dodgy" diffTableSolve, may also have run patterSolve or diffPat
  # If these worked, use them
  # If not, decide to contiue to fits or just use dodgy 
  
  if (usePattern & useDTP) 
  {
    # Favour DTP?
    usePattern = FALSE
    useTableVal3 = FALSE # Ditch dodgy
  } else if (usePattern) {
    useTableVal3 = FALSE # Ditch dodgy and use pattern
  } else if (useDTP) {
    useTableVal3 = FALSE # Ditch dodgy and use diff pattern
  } else if (useTableVal3 & acceptDodgy) {
    stop <- TRUE # Stop and leave useTableVal3 as TRUE
  } else if (useTableVal3 & !acceptDodgy) {
    useTableVal3 = FALSE # Ditch dodgy and continue to fits
    stop <- FALSE
  }
  
  # Try RR solve
  # Curentely - stops first score less than thresh
  # If using non-zero thresh, would be better to do all then find best
  if (!stop & RROn & m>2)
  {
    stop2 = FALSE
    p = 0
    while(!stop2)
    {
      p=p+1
      if (p==round((m-1)/2))
      {
        useRRVal <- FALSE
        stop2 = TRUE
      } else {
        res = RRSolve(sequ[1:(m-1)], p)
        sc = scorePCRR(res[[2]], sequ[m])
        # cat(paste("n:", n, "sc:", sc, "\n", sep=" "))
        if (!is.na(sc) & sc <= RRScore)
        {
          useRRVal <- TRUE
          stop <- TRUE
          stop2 <- TRUE
        }
      }
    }
  }
  
  if (!stop) # Do simple Fits
  {
    # Fit models using m points
    dFits = list()
    dScores = matrix(ncol=length(degs))
    dPreds = matrix(ncol=length(degs))
    dFits2 = list()
    dScores2 = matrix(ncol=length(degs))
    dPreds2 = matrix(ncol=length(degs))
    if (any(degs!=0))
    {
      for (d in degs)
      {
        dFits[[d]] <- fitModDeg(sequ, m-1, d, expFea=FALSE, plotOn)
        dScores[d] =  dFits[[d]][[2]]
        # dPreds[d] = dFits[[d]][[4]][m+1] 
        dFits2[[d]] <- fitModDeg(sequ, m-1, d, expFea=TRUE, plotOn)
        dScores2[d] =  dFits2[[d]][[2]]
        # dPreds2[d] = dFits2[[d]][[4]][m+1] 
      }
      # fits[[d]] = list with [[1]] = fit, [[2]] = score, [[3]] = yData, [[4]] = yPred
      
      # Find best score
      dBestScore = min(dScores)
      dBestScore2 = min(dScores2)
      dUseExp <- dBestScore2<dBestScore
      if (dUseExp)
      {
        dScores <- dScores2
        dBestScore <- min(dBestScore2)
      }
      dBestScore <- min(dBestScore)
    }
    
    # Fit models using m-1 points
    ppFits = list()
    ppScores = matrix(ncol=length(pps))
    ppFits2 <- ppFits
    ppScores2 <- ppScores
    if (any(pps!=0))
    {
      # Fit models using m-1 points
      ppFits = list()
      ppScores = matrix(ncol=length(pps))
      ppFits2 <- ppFits
      ppScores2 <- ppScores
      for (p in 1:length(pps))
      {
        if (m > 3*pps[p])
        {
          ppFits[[p]] <- fitModPP(sequ, m-1, pps[p], expFea=FALSE,  plotOn)
          ppScores[p] =  ppFits[[p]][[2]]
          ppFits2[[p]] <- fitModPP(sequ, m-1, pps[p], expFea=TRUE, plotOn)
          ppScores2[p] =  ppFits2[[p]][[2]]
          
        } else
        {
          ppScores[p] <- 9898
          ppScores2[p] <- 9898
        }
      }
      # fits[[d]] = list with [[1]] = fit, [[2]] = score, [[3]] = yData
      
      # Find best score and refit all data with that model
      ppBestScore = min(ppScores[!is.na(ppScores)])
      ppBestScore2 = min(ppScores2[!is.na(ppScores)])
      ppUseExp <- ppBestScore2<ppBestScore
      if (ppUseExp)
      {
        ppScores <- ppScores2
      }
      ppBestScore <- min(c(ppBestScore, ppBestScore2))
    }
    
    
    # Get the best fit score
    if (any(pps!=0) & any(degs!=0))
    {
      # Both kinds of fits done
      if (dBestScore < ppBestScore)
      {
        useFit <- 1
        useExp <- dUseExp
        bestScore = dBestScore
        useDegsFit <- TRUE
      }   else
      {
        useFit <- 2
        useExp <- ppUseExp
        bestScore = ppBestScore
        usePPSFit <- TRUE 
      }
      
    } else if (any(pps!=0)) # PPS on, degs off
    {
      useFit <- 2
      useExp <- ppUseExp
      usePPSFit <- TRUE
      bestScore = ppBestScore 
    } else if (any(degs!=0)) # Degs on, pps off
    {
      useFit <- 1
      useExp <- dUseExp
      bestScore = dBestScore
      useDegsFit <- TRUE
    }
    
    # Is this best score good enough?
    if ((bestScore < pcThresh) & !is.nan(bestScore))
    {
      stop = TRUE
    } else 
    {
      # Bad fit
      useDegsFit <- FALSE
      usePPSFit <- FALSE
      stop = FALSE
    }
  }
  
  svmScore = NaN
  # SVM1  
  if (!stop & degSVMOn)
  {
    # Table solve didn't work, simple fits didn't work. Try something non-parametric.
    
    svmFit <- tryCatch({
      svmFit = fitSVMDegs(sequ, m, plotOn=plotOn)
      
    }, error = function(err){
      svmFit = NaN 
      return(svmFit)
    }, finally = {
    })
    
    if (!length(svmFit)==1){
      svmScore = svmFit[[2]]
      if (svmScore<pcThreshSVM)
      {
        s <- "Used deg SVM"
        useSVM = TRUE
        stop = TRUE
        wwd = 7
      } 
    }
  }
  
  # SVM2  
  svmScore = NaN
  if (!stop & ppSVMOn)
  {
    
    svmFit <- tryCatch({
      pp = 1;
      svmFit = fitSVMPP(sequ, m, pp, plotOn = plotOn)
      
    }, error = function(err) {
      svmFit = NaN 
      return(svmFit)
    }, finally = {
      
    })
    
    # Table solve didn't work, simple fits didn't work. Try something non-parametric.
    if (!length(svmFit)==1) {
      svmScore = svmFit[[2]]
      if (svmScore<pcThreshSVM)
      {
        s <- "Used pp SVM"
        useSVM = TRUE
        stop = TRUE
        wwd = 8
      }
    }
  }
  
  ######################################################### What value to use?
  if (useRRVal)
  {
    # Predict next value, using whole sequence and p from above
    res = RRSolve(sequ, p)
    yPred[r,1] <- res[[2]]
    wwd = 12
    bestScore = sc
    s <- paste("Used RR val of degree", p)
  } else if (useDTP){
    resDTP = diffTablePattern(sequ)
    wwd = 16
    yPred[r,1] <- resDTP[[3]]
    yPred[r,2] <- wwd
    s <- paste("Used DTP with depth", resDTP[[5]][[1,1]], "and pattern at level", resDTP[[5]][[1,2]])
  } else if (useTableVal) {
    yPred[r,1] <- res[[2]][[1]]
    yPred[r,2] <- 3
    s <- "Used diff table"
  } else if (useTableVal2) {
    yPred[r,1] <- res[[2]][[1]]
    yPred[r,2] <- 30
    s <- paste(s, "Diff table solved with degree", dep, "(", res[[5]], ")", sep=" ")
  } else if (useTableVal3) {
    yPred[r,1] <- res[[2]][[1]]
    yPred[r,2] <- 31
    s <- paste(s, "Diff table solved with degree", dep, "(", res[[5]], ")", sep=" ")
  } else if (usePattern){
    yP = c(sequ, resPat[[2]])
    yPred[r,1] <- yP[m+1]
    yPred[r,2] <- wwd
    bestScore = resPat[[3]]
    s <- paste(s, "Pattern solved, score:", resPat[[3]], sep=" ")
  } else if (useDegsFit) {
    bestDeg <- degs[dScores==bestScore]
    
    if (length(bestDeg)>1)
    {
      # If multiple degrees score same, take lowest
      bestDeg <- bestDeg[1]
    }
    f <- fitModDeg(sequ, m, bestDeg, expFea=useExp, plotOn)
    
    # Predict next value
    newData <- f[[3]]
    newRow = nrow(newData)+1;
    newData[newRow,2:(bestDeg+1)] <- newRow^degs[1:bestDeg]
    
    
    s <- paste("used fit with degree", bestDeg, sep=" ")
    wwd = 3
    if (useExp)
    {
      newRowData = expandFeaturesNL(newData[newRow,1:(bestDeg+1)], "")
      newData[newRow,] <- newRowData[[1]]
      s <- paste(s, "+ extra features", sep=" ")
      wwd = 4
    } 
    
    # # Predict next point
    yPred[r,1] <- round(predict(f[[1]], newdata=data.frame(newData[newRow,])))
    yPred[r,2] <- 1
    
    # # Get already predicted next value NO
    # if (useExp)
    # {
    #   yPred[r,1] <- dFits2[[bestDeg]][[4]][length(dFits2[[bestDeg]][[4]])]
    #   s <- paste(s, "+ extra features", sep=" ")
    #   wwd = 3
    # } else {
    #   yPred[r,1] <- dFits[[bestDeg]][[4]][length(dFits[[bestDeg]][[4]])]
    #   wwd = 4
    # }
  } else if (usePPSFit) {
    bestpp <- pps[ppScores==bestScore]
    
    if (length(bestpp)>1)
    {
      # If multiple  score same, take lowest
      bestpp <- bestpp[1]
    }
    f <- fitModPP(sequ, m, bestpp, useExp, plotOn)
    
    # Predict next value
    newData <- f[[3]]
    newRow = nrow(newData)+1;
    # newData[newRow,2:(bestDeg+1)] <- newRow^degs[1:bestDeg]
    newData[newRow,1:(bestpp+1)] <- sequ[(m+1):((m+1)-(bestpp))]
    
    s <- paste("used fit with previous points", bestpp, sep=" ")
    wwd = 5
    if (useExp)
    {
      newRowData  = expandFeaturesNL(newData[newRow,1:(bestpp+1)], "")
      newData[newRow,] <- newRowData[[1]]
      s <- paste(s, "+ extra features", sep=" ")
      wwd = 6
    } 
    
    # Predict next point
    yPred[r,1] <- round(predict(f[[1]], newdata=data.frame(newData[newRow,])))
    yPred[r,2] <- 1
    # if (useExp)
    # {
    #   yPred[r,1] <- ppFits2[[bestpp]][[4]][length(ppFits2[[bestpp]][[4]])]
    #   s <- paste(s, "+ extra features", sep=" ")
    #   wwd = 5
    # } else {
    #   yPred[r,1] <- ppFits[[bestpp]][[4]][length(ppFits[[bestpp]][[4]])]
    #   wwd = 6
    # }
  } else if (useSVM) {
    yPred[r,1] <- svmFit[[4]][length(svmFit[[4]])]
    useFit = 3
    useExp = 0 # For now
    bestScore = svmFit[[2]]
  } else { # Fall back to mode
    # Save predicted value
    yPred[r,1] <- mode(sequ)
    # Save type of prediction for debugging
    # y Pred[r,2] <- 0
    s <- "Used mode"
    wwd=-1
  }
  
  # Record what was done
  yPred[r,4] <- useFit
  yPred[r,5] <- useExp
  yPred[r,3] <- bestScore
  yPred[r,2] <- wwd
  
  # Report
  pcd = paste("(", round(r/n*100, 2), "%)", sep="")
  cat("Done row", r, "/", n, pcd, s, "\n", sep=" ")
  
}

# Report number of sequences where a fit was used instead of the mode
sum(yPred[,2]==1)
# Report number of sequences where a table was used instead of the mode
sum(yPred[,2]==5)
# -1 (or NaN): Mode
sum(yPred[!is.na(yPred[,2]),2]==-1)
sum(is.na(yPred[,2]))
# 1212: RR
sum(yPred[!is.na(yPred[,2]),2]==12)
# 6565: DTP
sum(yPred[!is.na(yPred[,2]),2]==16)
# 1: diffTableSolve
sum(yPred[!is.na(yPred[,2]),2]==1)
# 2: diffTableSolve on seq sig
sum(yPred[!is.na(yPred[,2]),2]==2)
# 30: Difftable solve with extra degree
sum(yPred[!is.na(yPred[,2]),2]==30)
# 31: Difftable solve with extra degree (dodgy)
sum(yPred[!is.na(yPred[,2]),2]==31)
# 40: Pattern search
sum(yPred[!is.na(yPred[,2]),2]==40)
# 3: Deg fit
sum(yPred[!is.na(yPred[,2]),2]==3)
# 4: Deg fit with NL features
sum(yPred[!is.na(yPred[,2]),2]==4)
# 5: PP fit
sum(yPred[!is.na(yPred[,2]),2]==5)
# 6: PP fit with NL features
sum(yPred[!is.na(yPred[,2]),2]==6)
# 7: Deg SVM
sum(yPred[!is.na(yPred[,2]),2]==7)
# 8: PP (1) SVM
sum(yPred[!is.na(yPred[,2]),2]==8)


# Create data frame for submission
submission <- data.frame(Id=data$Id, Last=yPred[,1])

# Turn off scientific notation - does this really make a difference??
options(scipen=999)

# Write submission
write.csv(submission, "RSubmissionMaster19.csv", row.names=FALSE)