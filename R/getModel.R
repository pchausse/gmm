#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

getModel <- function(object, ...)
    {
        UseMethod("getModel")
    }

getModel.tsls <- function(object, ...)
    {
        obj <- getModel.baseGmm(object, ...)        
        return(obj)
    }

getModel.sysGmm <- function(object, ...)
    {
        if (object$commonCoef & !is.null(object$crossEquConst))
            {
                object$commonCoef <- FALSE
                warning("When crossEquConst is not NULL, commonCoef=TRUE is ignore and set to FALSE")
            }
        object$allArg <- c(object, ...)
        object$formula <- list(g=object$g, h=object$h)
        if (!is.list(object$g))
            stop("g must be list of formulas")
        if (length(object$g) == 1)
            stop("For single equation GMM, use the function gmm()")
        if (!all(sapply(1:length(object$g), function(i) is(object$g[[i]], "formula"))))
            stop("g must be a list of formulas")
        if (!is.list(object$h))
            {
                if(!is(object$h, "formula"))
                    stop("h is either a list of formulas or a formula")
                else
                    object$h <- list(object$h)
            } else {
                if (!all(sapply(1:length(object$h), function(i) is(object$h[[i]], "formula"))))
                    stop("h is either a list of formulas or a formula")
            }
        if (length(object$h) == 1)
            {
                object$h <- rep(object$h, length(object$g))
                typeDesc <- "System Gmm with common instruments"
                typeInst <- "Common"
            } else {
                if (length(object$h) != length(object$g))
                    stop("The number of formulas in h should be the same as the number of formulas in g, \nunless the instruments are the same in each equation, \nin which case the number of equations in h should be 1")                    
                typeDesc <- "System Gmm"
                typeInst <- "nonCommon"
            }
        if (object$commonCoef)
            typeDesc <- paste(typeDesc, " (Common Coefficients)")
        dat <- lapply(1:length(object$g), function(i) try(getDat(object$g[[i]], object$h[[i]], data = object$data,
                                                                 error=!object$commonCoef), silent=TRUE))
        chk <- sapply(1:length(dat), function(i) any(class(dat[[i]]) == "try-error"))
        chk <- which(chk)
        mess <- vector()
        for (i in chk)
            {
                mess <- paste(mess, "\nSystem of equations:", i, "\n###############\n", sep="")
                mess <- paste(mess, attr(dat[[i]], "condition")[[1]])
            }
        if (length(chk)>0)
            stop(mess)
        if (is.null(names(object$g)))
            names(dat) <- paste("System_", 1:length(dat), sep="")
        else
            names(dat) <- names(object$g)        
        object$gradv <- .DmomentFct_Sys
        object$formula <- list(g=object$g, h=object$h)
        if (!all(sapply(1:length(dat), function(i) dat[[i]]$ny == 1)))
            stop("The number of dependent variable must be one in every equation")
        
        clname <- "sysGmm.twoStep.formula"
        object$x <- dat
        namex <- lapply(1:length(dat), function(i) colnames(dat[[i]]$x[,2:(1+dat[[i]]$k), drop=FALSE]))
        nameh <- lapply(1:length(dat), function(i) colnames(dat[[i]]$x[,(2+dat[[i]]$k):(1+dat[[i]]$k+dat[[i]]$nh), drop=FALSE]))            
        namey <- lapply(1:length(dat), function(i) colnames(dat[[i]]$x[,1, drop=FALSE]))
        object$namesCoef <- namex
        namesgt <- lapply(1:length(dat), function(i) paste(namey[[i]], "_", nameh[[i]], sep=""))
        object$namesgt <- namesgt
        object$namesy <- namey
        attr(object$x,"ModelType") <- "linear"
       #for (i in 1:length(object$x))
        #    attr(object$x[[i]], c("linear")) <- attr(object$x, "modelType")
        attr(object$x, "k") <- lapply(1:length(dat), function(i) length(object$namesCoef[[i]]))
        attr(object$x, "q") <- lapply(1:length(dat), function(i) length(object$namesgt[[i]]))
        attr(object$x, "n") <- lapply(1:length(dat), function(i) nrow(object$x[[i]]$x))
        object$TypeGmm <- class(object)
        attr(object$x, "weight") <- list(w=object$weightsMatrix,
                                         centeredVcov=object$centeredVcov)
        attr(object$x, "weight")$WSpec <- list()
        attr(object$x, "weight")$WSpec$sandwich <- list(kernel = object$kernel, bw = object$bw,
                                                        prewhite = object$prewhite,
                                                        ar.method = object$ar.method,
                                                        approx = object$approx, tol = object$tol)
        attr(object$x, "weight")$vcov <- object$vcov
        k <- lapply(1:length(dat), function(i) dat[[i]]$k)
        q <- lapply(1:length(dat), function(i) dat[[i]]$nh)
        df <- lapply(1:length(dat), function(i) dat[[i]]$nh-dat[[i]]$k)
        k2 <- do.call(c,k)
        if (object$commonCoef | !is.null(object$crossEquConst))
            {
                if (!all(k2==k2[1]))
                    stop("In a common coefficient model the number of regressors is the same in each equation")
                if (object$commonCoef)
                    totK <- k2[1]
                else
                    totK <- length(dat)*(k2[1]-length(object$crossEquConst)) + length(object$crossEquConst)
                if (sum(do.call(c,q))<totK)
                    stop("The number of moment conditions is less than the number of coefficients")
                if (!is.null(object$crossEquConst))
                    {
                        object$crossEquConst <- sort(object$crossEquConst)
                        if (length(object$crossEquConst) == k2[1])
                            if (all(object$crossEquConst==(1:k2[1])))
                                {
                                    object$crossEquConst <- NULL
                                    object$commonCoef <- TRUE
                                }
                    }
            }
        attr(object$x, "sysInfo") <- list(k=k, df=df, q=q, typeInst=typeInst, typeDesc=typeDesc, commonCoef=object$commonCoef,
                                          crossEquConst=object$crossEquConst)
        object$g <- .momentFct_Sys
        class(object)  <- clname
        return(object)
    }


getModel.constGmm <- function(object, ...)
    {
        class(object) <- "baseGmm"
        obj <- getModel(object)
        if (!is.null(object$t0))
            {
                if (!is.null(dim(object$eqConst)))
                    stop("When t0 is provided, eqConst must be a vector which indicates which parameters to fix")
                if (length(object$eqConst)>=length(object$t0))
                    stop("Too many constraints; use evalGmm() if all coefficients are fixed")
                if (is.character(object$eqConst))
                    {
                        if (is.null(names(object$t0)))
                            stop("t0 must be a named vector if you want eqConst to be names")
                        if (any(!(object$eqConst %in% names(object$t0))))
                            stop("Wrong coefficient names in eqConst")
                        object$eqConst <- sort(match(object$eqConst,names(object$t0)))
                    }
                restTet <- object$t0[object$eqConst]
                obj$t0 <- object$t0[-object$eqConst]
                object$eqConst <- cbind(object$eqConst,restTet)
            }  else {
                if (is.null(dim(object$eqConst)))
                    stop("When t0 is not provided, eqConst must be a 2xq matrix")
            }
        attr(obj$x, "eqConst") <- list(eqConst = object$eqConst)
        rownames(attr(obj$x, "eqConst")$eqConst) <- obj$namesCoef[object$eqConst[,1]]
        object$eqConst <- attr(obj$x, "eqConst")$eqConst
        if(is(object$g, "formula"))
            {
                if (obj$x$ny>1)
                    stop("Constrained GMM not implemented yet for system of equations")
                if (obj$x$k<=0)
                    stop("Nothing to estimate")
            }
        obj$eqConst <- object$eqConst
        attr(obj$x, "k") <- attr(obj$x, "k")-nrow(object$eqConst)
        obj$namesCoef <- obj$namesCoef[-object$eqConst[,1]]
        obj$type <- paste(obj$type,"(with equality constraints)",sep=" ")	
        mess <- paste(rownames(object$eqConst), " = " , object$eqConst[,2], "\n",collapse="")
        mess <- paste("#### Equality constraints ####\n",mess,"##############################\n\n",sep="")
        obj$specMod <- mess
        return(obj)
    }

getModel.baseGmm <- function(object, ...)
    {
        object$allArg <- c(object, list(...))        
        if(is(object$g, "formula"))
            {
                object$gradv <- .DmomentFct
                object$gradvf <- FALSE
                dat <- getDat(object$g, object$x, data = object$data)
                if(is.null(object$weightsMatrix))
                    {
                        clname <- paste(class(object), ".", object$type, ".formula", sep = "")
                    } else {    
                        clname <- "fixedW.formula"
                        object$type <- "One step GMM with fixed W"
                    }
                object$x <- dat
                object$gform<-object$g
                namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k), drop=FALSE])
                nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh), drop=FALSE]) 
                if (dat$ny > 1)
                    {
                        namey <- colnames(dat$x[,1:dat$ny, drop=FALSE])
                        object$namesCoef <- paste(rep(namey, dat$k), "_",
                                                  rep(namex, rep(dat$ny, dat$k)), sep = "")
                        object$namesgt <- paste(rep(namey, dat$nh), "_",
                                                rep(nameh, rep(dat$ny, dat$nh)), sep = "")
                    } else {
                        object$namesCoef <- namex
                        object$namesgt <- nameh
                    }
                attr(object$x,"ModelType") <- "linear"
                attr(object$x, "k") <- object$x$k
                attr(object$x, "q") <- object$x$ny*object$x$nh
                attr(object$x, "n") <- NROW(object$x$x)
                attr(object$x, "namesgt") <- object$namesgt
            } else {
                attr(object$x,"ModelType") <- "nonlinear"
                attr(object$x, "momentfct") <- object$g
                if (object$optfct == "optimize")
                    attr(object$x, "k") <- 1
                else
                    attr(object$x, "k") <- length(object$t0)
                attr(object$x, "q") <- NCOL(gt <- object$g(object$t0, object$x))
                attr(object$x, "n") <- NROW(gt)
                if (object$optfct == "optimize")
                    {
                        object$namesCoef <- "Theta1"
                    } else {
                        if(is.null(names(object$t0)))
                            object$namesCoef <- paste("Theta[" ,1:attr(object$x, "k"), "]", sep = "")
                        else
                            object$namesCoef <- names(object$t0)
                    }
                if(is.null(object$weightsMatrix))
                    {
                        clname <- paste(class(object), "." ,object$type, sep = "")
                    } else {
                        clname <- "fixedW"
                        object$type <- "One step GMM with fixed W"
                        attr(object$x, "weight")$w <- object$weightsMatrix
                    }
                if (!is.function(object$gradv))
                    { 
                        object$gradvf <- FALSE
                    } else {
                        attr(object$x, "gradv") <- object$gradv    
                        object$gradvf <- TRUE
                    }
                object$gradv <- .DmomentFct
            }   
        object$TypeGmm <- class(object)
        attr(object$x, "weight") <- list(w=object$weightsMatrix,
                                         centeredVcov=object$centeredVcov)
        attr(object$x, "weight")$WSpec <- list()
        attr(object$x, "weight")$WSpec$sandwich <- list(kernel = object$kernel, bw = object$bw,
                                                        prewhite = object$prewhite,
                                                        ar.method = object$ar.method,
                                                        approx = object$approx, tol = object$tol)
        attr(object$x, "weight")$vcov <- object$vcov
        attr(object$x, "mustar") <- object$mustar
        object$g <- .momentFct
        class(object)  <- clname
        return(object)
    }

getModel.constGel <- function(object, ...)
    {
        class(object) <- "baseGel"
        obj <- getModel(object)
        if (!is.null(dim(object$eqConst)))
            stop("eqConst must be a vector which indicates which parameters to fix")
	if (length(object$eqConst)>=length(object$tet0))
            stop("Too many constraints; use evalGel() if all coefficients are fixed")
        if (is.character(object$eqConst))
            {
		if (is.null(names(object$tet0)))
                    stop("tet0 must be a named vector if you want eqConst to be names")
		if (any(!(object$eqConst %in% names(object$tet0))))
                    stop("Wrong coefficient names in eqConst")
		object$eqConst <- sort(match(object$eqConst,names(object$tet0)))
            }
	restTet <- object$tet0[object$eqConst]
	obj$tet0 <- object$tet0[-object$eqConst]
	object$eqConst <- cbind(object$eqConst,restTet)    
        attr(obj$x, "eqConst") <- list(eqConst = object$eqConst)
        rownames(attr(obj$x, "eqConst")$eqConst) <- obj$namesCoef[object$eqConst[,1]]
        object$eqConst <- attr(obj$x, "eqConst")$eqConst
        if(is(object$g, "formula"))
            {
                if (obj$x$ny>1)
                    stop("Constrained GMM not implemented yet for system of equations")
            }
        obj$eqConst <- object$eqConst
        attr(obj$x, "k") <- attr(obj$x, "k")-nrow(object$eqConst)
        obj$namesCoef <- obj$namesCoef[-object$eqConst[,1]]
        obj$typeDesc <- paste(obj$typeDesc,"(with equality constraints)",sep=" ")	
        mess <- paste(rownames(object$eqConst), " = " , object$eqConst[,2], "\n",collapse="")
        mess <- paste("#### Equality constraints ####\n",mess,"##############################\n\n",sep="")
        obj$specMod <- mess
        return(obj)
    }

getModel.baseGel <- function(object, ...)
    {
        object$allArg <- c(object, list(...))
        object$allArg$weights <- NULL
        object$allArg$call <- NULL
        if(is(object$g, "formula"))
            {
                dat <- getDat(object$g, object$x, data = object$data)
                k <- dat$k            
                if (is.null(object$tet0))
                    {
                        if (!is.null(object$eqConst))
                            stop("You have to provide tet0 with equality constrains")
                        if (object$optfct == "optimize")
                            stop("For optimize, you must provide the 2x1 vector tet0")
                        res0 <- gmm(object$g, object$x, data=object$data)
                        object$tet0 <- res0$coefficients
                        if (object$smooth)
                            gt <- res0$gt
                    } else {
                        if (object$optfct == "optimize")
                            {
                                if (k != 1)
                                    stop("optimize() is for univariate optimization")
                                if (length(object$tet0) != 2)
                                    stop("For optimize(), tet0 must be a 2x1 vector")
                            } else {
                                if (k != length(object$tet0))
                                    stop("The number of starting values does not correspond to the number of regressors")
                            }                
                        if (object$smooth)
                            gt <- gmm(object$g, object$x, data=object$data)$gt
                    }                
                clname <- paste(class(object), ".modFormula", sep = "")
                object$gradv <- .DmomentFct
                object$gradvf <- FALSE
                object$x <- dat
                object$gform<-object$g
                namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k), drop=FALSE])
                nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh), drop=FALSE])
                if (dat$ny > 1)
                    {
                        namey <- colnames(dat$x[,1:dat$ny])
                        namesCoef <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
                        object$namesgt <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
                    } else {
                        namesCoef <- namex
                        object$namesgt <- nameh
                    }
                if (is.null(names(object$tet0)))
                    object$namesCoef <- namesCoef
                else
                    object$namesCoef <- names(object$tet0)
                attr(object$x,"ModelType") <- "linear"
                attr(object$x, "k") <- k
                attr(object$x, "q") <- object$x$ny*object$x$nh
                
                attr(object$x, "n") <- NROW(object$x$x)
            } else {
                if (is.null(object$tet0))
                    stop("You must provide the starting values tet0 for nonlinear moments")
                if(any(object$optfct == c("optim", "nlminb")))
                    k <- length(object$tet0)
                else
                    k <- 1                    
                attr(object$x,"ModelType") <- "nonlinear"
                attr(object$x, "momentfct") <- object$g
                attr(object$x, "k") <- k
                attr(object$x, "q") <- NCOL(gt <- object$g(object$tet0, object$x))
                attr(object$x, "n") <- NROW(gt)
                if(is.null(names(object$tet0)))
                    object$namesCoef <- paste("Theta[" ,1:attr(object$x, "k"), "]", sep = "")
                else
                    object$namesCoef <- names(object$tet0)
                if (!is.function(object$gradv) | object$smooth)
                    { 
                        object$gradvf <- FALSE
                    } else {
                        attr(object$x, "gradv") <- object$gradv    
                        object$gradvf <- TRUE                            
                    }
                object$gradv <- .DmomentFct
                if (object$smooth)
                    gt <- gmm(object$g, object$x, object$tet0, wmatrix = "ident", ...)$gt
                clname <- paste(class(object), ".mod", sep = "")
            }
        if (object$smooth)
            {
                if (is.function(object$gradv))
                    warning("The provided gradv is not used when smooth=TRUE",
                            call. = FALSE)		
                if(object$kernel == "Truncated")
                    {
                        object$wkernel <- "Bartlett"
                        object$k1 <- 2
                        object$k2 <- 2
                    }
                if(object$kernel == "Bartlett")
                    {
                        object$wkernel <- "Parzen"
                        object$k1 <- 1
                        object$k2 <- 2/3
                    }
                gt <- scale(gt, scale=FALSE)
                class(gt) <- "gmmFct"
                if (is.function(object$bw))
                    object$bwVal <- object$bw(gt, kernel = object$wkernel, prewhite = object$prewhite, 
                                              ar.method = object$ar.method, approx = object$approx)
                else
                    object$bwVal <- object$bw
                object$w <- smoothG(gt, bw = object$bwVal, kernel = object$kernel, tol = object$tol_weights)$kern_weights
                attr(object$x,"smooth") <- list(bw=object$bwVal, w=object$w, kernel=object$kernel)
            } else {
                object$k1 <- 1
                object$k2 <- 1
                object$w <- kernel(1)
                object$bwVal <- 1
            }
        object$g <- .momentFct
        object$CGEL <- object$alpha
        object$typeDesc <- object$type
        class(object) <- clname
        return(object)
    }

getModel.ateGel <- function(object, ...)
    {
        object$allArg <- c(object, list(...))
        object$allArg$weights <- NULL
        object$allArg$call <- NULL
        if(is(object$g, "formula"))
            {
                dat <- getDat(object$g, object$x, data = object$data)
                if (!is.null(object$w))
                    if (is(object$w, "formula"))
                        {
                            dat$w <- model.matrix(object$w, object$data)[,-1,drop=FALSE]
                        } else {
                            stop("w must be a formula")
                        }                
                if (dat$ny>1 | dat$ny==0)
                    stop("You need one and only one dependent variable")
                k <- dat$k
                if (k>2 & object$momType=="ATT")
                    stop("Cannot compute ATT with multiple treatments")
                if (attr(dat$mt, "intercept")!=1)
                    stop("An intercept is needed to compute the treatment effect")
                if (!all(dat$x[,3:(k+1)] %in% c(0,1)))
                    stop("The treatment indicators can only take values 0 or 1")
                if (colnames(dat$x)[k+2] == "(Intercept)")
                    {
                        dat$x <- dat$x[,-(k+2)]
                        dat$nh <- dat$nh-1
                    }
                if (!is.null(object$popMom))
                    {
                        if (length(object$popMom)!=dat$nh)
                            stop("The dim. of the population moments does not match the dim. of the vector of covariates")                                 
                    }
                if (is.null(object$tet0))
                    {
                        if (is.null(dat$w))
                            {
                                tet0 <- lm(dat$x[,1]~dat$x[,3:(k+1)])$coef
                            } else {
                                tet0 <- lm(dat$x[,1]~dat$x[,3:(k+1)]+dat$w)$coef
                            }
                        tet0 <- c(tet0, colMeans(dat$x[,3:(k+1),drop=FALSE]))
                        names(tet0) <- NULL
                        object$tet0 <- tet0
                    } else {
                        ntet0 <- 2*k-1 + ifelse(is.null(dat$w), 0, ncol(dat$w))
                        if (length(object$tet0) != ntet0)
                            stop("tet0 should have a length equal to 2x(number of treatments)+1+number of w's if any")
                    }
                if (object$family != "linear")
                    {
                        if (any(!(dat$x[,1]%in%c(0,1))))
                            stop("For logit or probit, Y can only take 0s and 1s")
                        family <- binomial(link=object$family)
                        if (object$family == "logit")
                            family$mu.eta2 <- function(x, family) family$mu.eta(x)*(1-2*family$linkinv(x))
                        else
                            family$mu.eta2 <- function(x, family) -x*family$mu.eta(x)
                        
                    } else {
                        family <- NULL
                    }
                q <- dat$nh + 2*k+1
                if (object$momType != "bal" | !is.null(object$popMom))
                    q  <- q+dat$nh
                if (!is.null(dat$w))
                    {
                        q <- q+ncol(dat$w)
                        namew <- colnames(dat$w)
                    } else {
                        namew <- NULL
                    }
                object$gradv <- .DmomentFctATE
                object$x <- dat
                object$gradvf <- FALSE
                object$gform<-object$g                
                namex <- colnames(dat$x[,2:(k+1)])
                nameh <- colnames(dat$x[,(k+2):ncol(dat$x), drop=FALSE])
                namesCoef <- c(namex, namew, paste("TreatProb", 1:(k-1), sep=""))
                namesgt <- paste(rep(paste("Treat", 1:(k-1), sep=""),
                                     rep(dat$nh, k-1)), "_", rep(nameh, k-1), sep="")
                object$namesgt <- c(namesCoef,namesgt)
                if (object$momType != "bal" | !is.null(object$popMom))
                    object$namesgt <- c(object$namesgt, paste("bal_", nameh, sep=""))
                if (is.null(names(object$tet0)))
                    object$namesCoef <- namesCoef
                else
                    object$namesCoef <- names(object$tet0)
                attr(object$x,"ModelType") <- "linear"
                attr(object$x, "k") <- k
                attr(object$x, "q") <- q
                attr(object$x, "n") <- nrow(dat$x)
                attr(object$x, "momType") <- object$momType
                attr(object$x, "popMom") <- object$popMom
                attr(object$x, "family") <- family
            } else {
                stop("Not implemented yet for nonlinear regression")
            }
        if (object$momType == "ATT")
            metname <- "ATT"
        else
            metname <- "ATE"
        if (!is.null(object$popMom))
            {
                metname2 <- " with balancing based on population moments"
            } else {
                if (object$momType == "balSample")
                    metname2 <- " with balancing based on the sample moments"
                else
                    metname2 <- " with unrestricted balancing"
            }
        metname3 <- paste("\nMethod: ", object$type, sep="")
        if (!is.null(family))
            metname3 <- paste(metname3, ", Family: Binomial with ", family$link, " link", sep="")
        clname <- "baseGel.mod"
        object$k1 <- 1
        object$k2 <- 1
        object$w <- kernel(1)
        object$bwVal <- 1
        object$g <- .momentFctATE
        object$CGEL <- object$alpha
        object$typeDesc <- paste(metname, metname2, metname3, sep="")
        class(object) <- clname
        return(object)
    }

