library(rpart)
#source("honestRpart.R")

OUT_LOG_FLAG <- TRUE
OUT_LOG_LEVEL <- 0
GLOBAL_MAX <- 1e100

logging <- function(x=NULL, msg="", level=1) {
    if(OUT_LOG_FLAG && level > OUT_LOG_LEVEL) {
        cat(paste("Logging << ", msg, "\n"))
        if(!is.null(x)) {
            print(x)
        }
    }
}

colorize <- function(x) {
    x <- rank(x)
    if (max(x) == min(x)) {
        return(rep(0.5, length(x)))
    } else {
        return((x-min(x)) / (max(x)-min(x)))
    }
}
#### Witchwoods Part ####
crown <- function(z) {
    colnames(z) <- paste("z", 1:ncol(z), sep="")
    return(z)
}

pred.woods <- function(z, model) {
    if (is.vector(z)) {
        z <- matrix(z, nrow=1)
    }
    tmpBeta <- matrix(0, nrow=nrow(z), ncol=model$p)
    for (i in 1:p) {
        tmpBeta[, i] <- predict(model$woods[[1]][[i]], crown(z))
        if (model$woods.len > 1) {
            for (j in 2:model$woods.len) {
                tmpBeta[, i] <- tmpBeta[, i] + model$lambda *  predict(model$woods[[j]][[i]], crown(z))
            }
        }
    }
    return(tmpBeta)
}

#### for ols
ols.diff <- function(y, x, beta) {
    return((y - rowSums(x * beta)) * x)
}

ols.init.beta <- function(y, x) {
    return(matrix(as.vector(lm(y ~ x+0)$coef), 
                  nrow=nrow(x), ncol=ncol(x), byrow=TRUE))
}

ols.pred <- function(x, beta) {
    return(rowSums(x * beta))
}

ols.loss <- function(y1, y2) {
    return(mean((y1-y2)^2))
}

#### for logistic
sigmoid <- function(z) {
    return(1/(1+exp(-z)))    
}

lg.pred <- function(x, beta) {
    return(rowSums(x * beta))
}

lg.diff <- function(y, x, beta) {
    return(x*(y-sigmoid(lg.pred(x, beta))))
}

lg.init.beta <- function(y, x) {
    return(matrix(as.vector(glm(y ~ x+0,
                                family = binomial(link = "logit"))$coef), 
                  nrow=nrow(x), ncol=ncol(x), byrow=TRUE))
}

lg.loss <- function(y1, y2) {
    return(-mean((1-y1)*(-y2)-log(1+exp(-y2))))
}

grow.woods <- function(u, z, control, is.plinear = FALSE) {
    if (ncol(u) != model$p) {
        logging(msg="Err: dim mismatch")
    }
    if (is.plinear) {
        w_control <- control
        w_control$cp <- 1e10
    }
    tmpWoods <- list()
    if (is.plinear) {
        tmpWoods[[1]] <- rpart(u[, 1] ~ crown(z), control=control, y = FALSE)
        if (model$p > 1) {
            for (i in 2:model$p) {
                tmpWoods[[i]] <- rpart(u[, i] ~ crown(z), control=w_control, y = FALSE)
            }
        }
    } else {
        for (i in 1:model$p) {
            tmpWoods[[i]] <- rpart(u[, i] ~ crown(z), control=control, y = FALSE)
        }
    }
    return(tmpWoods)
}

grow.woods.v2 <- function(u, z, model) {
    if (ncol(u) != model$p) {
        logging(msg="Err: dim mismatch")
    }
    if (model$method == "plinear") {
        w_control <- model$control
        w_control$cp <- 1e10
    }
    tmpSubset <- sample(model$n, model$n * model$subsample)
    tmpDF <- as.data.frame(crown(z))
    tmpWoods <- list()
    if (model$method == "plinear") {
        tmpDF$response = u[, 1]
        tmpWoods[[1]] <- rpart(response ~ ., data=tmpDF, subset=tmpSubset, control=model$control, y = FALSE)
        if (model$p > 1) {
            for (i in 2:model$p) {
                tmpDF$response = u[, i]
                tmpWoods[[i]] <- rpart(response ~ ., data=tmpDF, subset = tmpSubset, control=w_control, y = FALSE)
            }
        }
    } else if (model$method == "boulevard") {
        for (i in 1:model$p) {
            tmpWoods[[i]] <- honest.rpart(Y=u[, i], X=z, structY=u[, i], 
                                          subset=tmpSubset, 
                                          control=model$control)
        }
    } else {
        for (i in 1:model$p) {
            tmpDF$response = u[, i]
            tmpWoods[[i]] <- rpart(response ~., data=tmpDF, subset=tmpSubset, control=model$control, y = FALSE)
        }
    }
    return(tmpWoods)
}

predict.woods <- function(tree, model, newdata = NULL, train=FALSE) {
    if (model$method == "boulevard") {
        return(honest.rpart.predict(tree, newdata = newdata))
    # } else if (train) {
    #     return(predict(tree))
    } else {
        return(predict(tree, newdata = as.data.frame(crown(newdata))))
    }
}

lgd <- function(y, x, z, model) {
    #### set default values
    if (is.null(model$savetrees)) {
        model$savetrees = TRUE
    }
    if (is.null(model$xscale)) {
        model$xscale = FALSE
    }
    if (is.null(model$yscale)) {
        model$yscale = FALSE
    }
    if (is.null(model$dummy)) {
        model$dummy = FALSE
    }
    if (is.null(model$method)) {
        model$method = "ordinary"
    } 
    if (is.null(model$quantile_cut)) {
        model$quantile_cut = 0.05
    }
    if (model$method == "boulevard") {
        if (is.null(model$control$minsplits)) {
            model$control$minsplits = 10
        }
        if (is.null(model$subsample)) {
            model$subsample = 0.5
        }
    }
    
    #### scale predictor / responses
    if (model$xscale) {
        x <- scale(x)
        # print(head(x))
        model$xscale.center = attr(x, "scaled:center")
        model$xscale.scale = attr(x, "scaled:scale")
        if (model$dummy) {
            x <- cbind(1, x)
            model$xscale.center <- c(-1, model$xscale.center)
            model$xscale.scale <- c(1, model$xscale.scale)
        }
        if (model$yscale) {
            y <- scale(y)
            model$yscale.center = attr(y, "scaled:center")
            model$yscale.scale = attr(y, "scaled:scale")
            y <- as.vector(y)
        } else {
            model$yscale.center = 0
            model$yscale.scale = 1
        }
    } else if (model$dummy) {
        x <- cbind(1, x)
    }
    
    if (model$method == "boulevard") {
        x <- x * sqrt(model$n)
    }
    z <- crown(z)
    #### lgb
    model$lc <- c()
    # if (model$method == "boulevard") {
    #     tmpBeta <- model$diff(y, x, 0)
    # } else {
        tmpBeta <- model$init(y, x)
    # }
    model$woods[[1]] <- grow.woods.v2(tmpBeta, z, model)
    model$beta <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    for (i in 1:model$p) {
        model$beta[, i] <- predict.woods(model$woods[[1]][[i]], model, newdata=z)
    }
    model$yhat <- model$pred(x, model$beta)
    model$lc <- c(model$lc, model$loss(y, model$yhat))
    if (model$ntree > 1) {
        for (ntree in 2:model$ntree) {
            gc()
            if (model$savetrees) {
                w_ntree = ntree
            } else {
                w_ntree = 1
            }
            #if (ntree %% 5 == 0) {cat("ITER:", ntree, "\tSave to", w_ntree)} #[modification by SKD to suppress print statements]
            # if (model$method == "boulevard") {
            #     u <- model$diff(y, x, model$beta / ntree * (ntree - 1))
            # } else {
                u <- model$diff(y, x, model$beta)
            # }
            for (i in 1:model$p) {
                if (!is.null(model$gradient_truncate)) {
                    quantile_truncate <- c(-model$gradient_truncate, model$gradient_truncate)
                } else {
                    quantile_truncate = as.vector(quantile(u[, i], 
                                                           probs=c(model$quantile_cut, 1-model$quantile_cut)))
                    u[u[, i] < quantile_truncate[1], i] = quantile_truncate[1]
                    u[u[, i] > quantile_truncate[2], i] = quantile_truncate[2]
                }
            }
            model$woods[[w_ntree]] <- grow.woods.v2(u, z, model)
            for (i in 1:model$p) {
                if (model$method == "ordinary") {
                    model$beta[, i] <- model$beta[, i] + model$lambda * predict.woods(model$woods[[w_ntree]][[i]], model, newdata=z)
                } else if (model$method == "boulevard") {
                    model$beta[, i] <- model$beta[, i] * (ntree - 1) / ntree + 
                        model$lambda / ntree * predict.woods(model$woods[[w_ntree]][[i]], model, newdata=z)
                } else {
                    model$beta[, i] <- model$beta[, i] + model$lambda * predict.woods(model$woods[[w_ntree]][[i]], model, newdata=z)
                }
            }
            # if (ntree %% 10 == 0) {print(head(model$beta))}
            model$yhat <- model$pred(x, model$beta)
            # plot(y, model$yhat)
            if (model$method=="boulevard") {
                tmpInflate <- lm(y ~ model$pred(x, model$beta) + 0)$coef[1]
                # boxplot(y / model$pred(x, model$beta), breaks=50)
                #if(ntree %% 5 == 0) {cat("\tblvRefLOSS:", model$loss(y, model$pred(x, model$beta * (1 + model$lambda) / model$lambda)))}
                #if(ntree %% 5 == 0) {cat("\tInflate:", tmpInflate)}
                #if(ntree %% 5 == 0) {cat("\tblvInfLOSS:", model$loss(y, model$pred(x, model$beta * tmpInflate)))}
                tmpLoss <- model$loss(y, model$yhat)
            } else {
                tmpLoss <- model$loss(y, model$yhat)
            }
            # tmpLoss <- model$loss(y, model$yhat)
            model$lc <- c(model$lc, tmpLoss)
            #if(ntree %% 5 == 0) {cat("\tLOSS:", tmpLoss, "\n")}
        }
    }
    
    if (model$method == "boulevard") {
        if (!is.null(model$inflate)) {
            model$inflate = lm(y ~ model$pred(x, model$beta) + 0)$coef[1]
        }
        model$beta <- model$beta * sqrt(model$n)
    }
    #### Unscale
    if (model$xscale) {
        for (i in 1:model$p) {
            model$beta[, i] <- model$beta[, i] * model$yscale.scale / model$xscale.scale[i]
        }
        if (model$dummy) model$beta[, 1] <- -model$beta%*%model$xscale.center  + model$yscale.center
    }
    return(model)
}

lgd.predict <- function(x, z, model, flag = FALSE) {
    if (model$dummy) {
        x <- cbind(1, x)
    }
    # print(head(x))
    tmpBeta <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    for (ntree in 1:model$ntree) {
        for (i in 1:model$p) {
            if (model$method == "ordinary") {
                if (ntree == 1) {
                    tmpBeta[, i] <- predict.woods(model$woods[[ntree]][[i]], model, newdata = crown(z))
                } else {
                    tmpBeta[, i] <- tmpBeta[, i] + model$lambda * predict.woods(model$woods[[ntree]][[i]], model, newdata = crown(z))
                }
            } else if (model$method == "boulevard") {
                tmpBeta[, i] <- tmpBeta[, i] * (ntree-1) / ntree + model$lambda / ntree * predict.woods(model$woods[[ntree]][[i]], model, newdata = crown(z))
            }
        }
    }
    if (model$method == "boulevard") {
        tmpBeta <- tmpBeta * sqrt(model$n)
    }
    #### Unscale
    if (model$xscale) {
        for (i in 1:model$p) {
            tmpBeta[, i] <- tmpBeta[, i] * model$yscale.scale / model$xscale.scale[i]
        }
        if (model$dummy) tmpBeta[, 1] <- -tmpBeta%*%model$xscale.center  + model$yscale.center
    }
    if (model$method == "boulevard") {
        if (flag) {
            tmpBeta = tmpBeta * (1 + model$lambda) / model$lambda 
        }
    }
    return(model$pred(x, tmpBeta))
}
