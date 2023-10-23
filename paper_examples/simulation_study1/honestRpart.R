library(rpart)
library(treeClust)
#library(plotrix)

honest.rpart.structure <- function(X, Y, method="standard", structY=NULL, leaf.size=3, control=NULL, diameter.test=NULL) {
  n <- nrow(X)
  colnames(X) <- paste("x", 1:ncol(X), sep="")
  pass.diameter.check <- FALSE
  while (!pass.diameter.check) {
    if (method=="random") {
      trainData <- data.frame(Y=rnorm(n), X=X)
    } else if (method == "hybrid") {
      trainData <- data.frame(Y = rnorm(n, 0, sd(structY)) + structY, X = X)
    } else {
      trainData <- data.frame(Y=structY, X=X)
    }
    if (is.null(control)) {
        tree.structure <- rpart(Y ~ . , data=trainData, control=rpart.control(cp = 0, minsplit=leaf.size+1), y = FALSE)
    } else {
        tree.structure <- rpart(Y ~ . , data=trainData, control=control, y = FALSE)
    }
    if (!is.null(diameter.test)) {
      colnames(diameter.test) <- paste("x", 1:ncol(X), sep="")
      where <- rpart.predict.leaves(tree.structure, newdata=data.frame(newdata), type="where")
      if (length(where)==length(unique(where))) return(tree.structure)
    } else {
      return(tree.structure)
    }
    method <- "random"
  }
}

honest.rpart <- function(X, Y, method="standard", structY=NULL, subset=NULL, leaf.size=3, control = NULL, diameter.test=NULL) {
  colnames(X) <- paste("x", 1:ncol(X), sep="")
  trainData <- data.frame(Y, X)
  
  tree <- list()
  tree$train.where <- rep(0, length(Y))
  if (is.null(subset)) {
    newX <- X
    newY <- Y
  } else {
    subset <- sort(subset)
    newX <- X[subset, ]
    newY <- Y[subset]
  }
  if (!is.null(diameter.test)) {
    tree$rpart.tree <- honest.rpart.structure(X, Y, method=method, structY=structY, leaf.size=leaf.size, control=control, diameter.test=diameter.test)
  } else {
    tree$rpart.tree <- honest.rpart.structure(X, Y, method=method, structY=structY, leaf.size=leaf.size, control=control)
  }
  
  tmpWhere <- unique(tree$rpart.tree$where)
  where <- rpart.predict.leaves(tree$rpart.tree, newdata=data.frame(X=newX), type="where")
  tmpPredict <- c()
  for (i in tmpWhere) {
    tmpCount <- sum(where==i)
    if (tmpCount > 0) {
      tmpPredict <- c(tmpPredict, mean(newY[where==i]))
      tree$train.where[subset[where==i]] = i
    } else {
      tmpPredict <- c(tmpPredict, 0)
    }
  }
  tree$predict <- tmpPredict
  names(tree$predict) <- tmpWhere
  return(tree)
}

honest.rpart.predict <- function(tree, newdata) {
  if (is.vector(newdata)) {
    newdata = matrix(newdata, nrow=1)
  }
  colnames(newdata) = paste("x", 1:ncol(newdata), sep="")
  where <- rpart.predict.leaves(tree$rpart.tree, newdata=data.frame(X=newdata), type="where")
  return(as.vector(tree$predict[paste(where)]))
}

honest.rpart.predict.weight <- function(tree, newdata) {
  if (is.vector(newdata)) {
    newdata = matrix(newdata, nrow=1)
  }
  colnames(newdata) = paste("x", 1:ncol(newdata), sep="")
  where <- rpart.predict.leaves(tree$rpart.tree, newdata=data.frame(X=newdata), type="where")
  ans <- t(sapply(where, function(x) return(tree$train.where==x)))*1
  for (i in 1:nrow(ans)){
    w <- sum(ans[i, ])
    if(w>0) {
      ans[i, ] = ans[i, ] / w
    }
  }
  return(ans)
}

predict.boulevard <- function(blv, X) {
  ntree <- length(blv$trees)
  lambda <- blv$lambda
  ans <- rep(0, nrow(X))
  for (b in 1:ntree) {
    ans <- (b-1)/b*ans + lambda/b*honest.rpart.predict(blv$trees[[b]], newdata = X)
  }
  return(ans*(1+lambda)/lambda)
}

predict.boulevard.variance <- function(blv, newdata, narrow = FALSE) {
  ntree <- length(blv$trees)
  lambda <- blv$lambda
  ans <- honest.rpart.predict.weight(blv$trees[[1]], newdata=newdata)
  if (ntree > 1) {
    for (b in 2:ntree) {
      ans <- ans + honest.rpart.predict.weight(blv$trees[[b]], newdata=newdata)
    }
    ans <- ans / ntree
  }
  ans <- apply(ans, 1, function(x) sum(x*x))
  ans <- ans * (1+lambda)^2 * blv$sigma2
  if (narrow) {
    ans <- ans / (1+lambda)^2
  } 
  return(ans)
}

boulevard <- function(X, Y, ntree=1000, lambda = 0.8, subsample=0.8, xtest=NULL, ytest=NULL, leaf.size=10, method="random") {
  n <- nrow(X)
  tree <- list()
  nss <- floor(n * subsample)
  ans <- rep(0, nrow(X))
  if (!is.null(xtest)) {
    predtest <- rep(0, nrow(xtest))
  }
  trainmse <- c()
  testmse <- c()
  for (b in 1:ntree) {
    if (b%%50 == 0) {
      cat("Training Iteration:", b, "\n")
    }
    res <- Y - ans
    if (method=="random") {
      tree[[b]] <- honest.rpart(X, res, method=method, subset=sample(n, nss, replace=FALSE), leaf.size=leaf.size)
    } else if (method=="standard") {
      tree[[b]] <- honest.rpart(X, res, method=method, structY=res, subset=sample(n, nss, replace=FALSE), leaf.size=leaf.size)
    }
    ans <- (b-1)/b*ans + lambda/b*honest.rpart.predict(tree[[b]], newdata = X)
    trainmse <- c(trainmse,mean((ans/lambda*(1+lambda)-Y)^2))
    if (!is.null(xtest)) {
      predtest <- (b-1)/b*predtest + lambda/b*honest.rpart.predict(tree[[b]], newdata=xtest)
      testmse <- c(testmse, mean((predtest/lambda*(1+lambda)-ytest)^2))
    }
  }
  return(list(trees=tree, 
              mse=trainmse, 
              testmse=testmse, 
              lambda=lambda,
              sigma2=trainmse[ntree]))
}


