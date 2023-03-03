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


.invTest <- function(object, which, level = 0.95, fact = 3, type=c("LR", "LM", "J"), corr=NULL)
    {
        type <- match.arg(type)
        argCall <- object$allArg
        if (length(which) > 1)
            stop("tests are inverted only for one parameter")        
        if (is.character(which))
            {
             which <- which(which == names(object$coefficients))   
             if (length(which) == 0)
                 stop("Wrong parameter names")
            }
        wtest <- which(type==c("LR", "LM", "J"))
        if (object$df > 0)
            {
                test0 <- specTest(object)$test[wtest,]
                test0 <- test0[1]
            } else {
                test0 <- 0
            }
        if (length(object$coefficients) == 1)
            {
                argCall$optfct <- NULL
                argCall$constraint <- NULL
                argCall$eqConst <- NULL
                argCall$eqConstFullVcov <- NULL
                fctCall <- "evalGel"
            } else {
                fctCall <- "gel"
                argCall$eqConst <- which
            }
        if (length(object$coefficients) == 2)
            {
                sdcoef <- sqrt(diag(vcov(object))[-which])
                coef <- object$coefficients[-which]
                upper <- coef + fact*sdcoef
                lower <- coef - fact*sdcoef
                argCall$method <- "Brent"
                argCall$lower <- lower
                argCall$upper <- upper
            }
        sdcoef <- sqrt(diag(vcov(object))[which])
        coef <- object$coefficients[which]
        int1 <- c(coef, coef + fact*sdcoef)
        int2 <- c(coef - fact*sdcoef, coef)
        fct <- function(coef, which, wtest, level, test0, corr=NULL)
            {
                argCall$tet0 <- object$coefficients
                argCall$tet0[which] <- coef
                obj <- do.call(get(fctCall), argCall)
                test <- as.numeric(specTest(obj)$test[wtest,1]) - test0
                if (is.null(corr))
                    level - pchisq(test, 1)
                else
                   level - pchisq(test/corr, 1)
            }
        res1 <- uniroot(fct, int1, which = which, wtest=wtest, level=level,
                        test0=test0, corr=corr)
        res2 <- uniroot(fct, int2, which = which, wtest=wtest, level=level,
                        test0=test0, corr=corr)
        sort(c(res1$root, res2$root))
    }
        

confint.gel <- function(object, parm, level = 0.95, lambda = FALSE,
                        type = c("Wald", "invLR", "invLM", "invJ"),
                        fact = 3, corr = NULL, ...)
    {
        type <- match.arg(type)
        z <- object	
        n <- nrow(z$gt)
        if (type == "Wald")
            {
                ntest <- "Direct Wald type confidence interval"
                se_par <- sqrt(diag(z$vcov_par))
                par <- z$coefficients
                tval <- par/se_par
                se_parl <- sqrt(diag(z$vcov_lambda))
                lamb <- z$lambda                
                zs <- qnorm((1 - level)/2, lower.tail=FALSE)
                ch <- zs*se_par                
                if(!lambda)
                    {
                        ans <- cbind(par-ch, par+ch)
                        dimnames(ans) <- list(names(par), c((1 - level)/2, 0.5+level/2))
                    }
                if(lambda)
                    {
                        if (length(z$coefficients) == length(z$lambda))
                            {
                                cat("\nNo confidence intervals for lambda when the model is just identified.\n")
                                return(NULL)
                            } else {
                                chl <- zs*se_parl
                                ans <- cbind(lamb - chl, lamb + chl)
                                dimnames(ans) <- list(names(lamb), c((1 - level)/2, 0.5 + level/2))
                            }
                    }		
                if(!missing(parm))
                    ans <- ans[parm,]
            } else {
                if(missing(parm))
                    parm <- names(object$coefficients)
                type <- strsplit(type, "v")[[1]][2]
                ntest <- paste("Confidence interval based on the inversion of the ", type, " test", sep="")
                ans <- lapply(parm, function(w) .invTest(object, w, level = level, fact = fact, type=type, corr=corr))
                ans <- do.call(rbind, ans)
                if (!is.character(parm))
                    parm <- names(object$coefficients)[parm]
                dimnames(ans) <- list(parm, c((1 - level)/2, 0.5+level/2))
            }    
        ans <- list(test=ans,ntest=ntest)
        class(ans) <- "confint"
        ans
    }

print.confint <- function(x, digits = 5, ...)
    {
        cat("\n", x$ntest, sep="")
        cat("\n#######################################\n")
	print.default(format(x$test, digits = digits),
                      print.gap = 2, quote = FALSE)
        invisible(x)
    }

coef.gel <- function(object, lambda = FALSE, ...) 
    {
	if(!lambda)
            object$coefficients
	else
            object$lambda
    }

vcov.gel <- function(object, lambda = FALSE, ...) 
    {
	if(!lambda)
            object$vcov_par
	else
            object$vcov_lambda
    }

print.gel <- function(x, digits = 5, ...)
    {
	if (is.null(x$CGEL))
            cat("Type de GEL: ", x$typeDesc, "\n")
	else
            cat("CGEL of type: ", x$typeDesc, " (alpha = ", x$CGEL, ")\n")
	if (!is.null(attr(x$dat,"smooth")))
            {
		cat("Kernel: ", attr(x$dat,"smooth")$kernel," (bw=",
                    attr(x$dat,"smooth")$bw,")\n\n")
            }
	else
            cat("\n")
        
	cat("Coefficients:\n")
	print.default(format(coef(x), digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\n")
        cat("Lambdas:\n")
        print.default(format(coef(x, lambda = TRUE), digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\n")
	cat("Convergence code for the coefficients: ", x$conv_par,"\n")
        cat("Convergence code for Lambda: ", x$conv_lambda$convergence,"\n")
        cat(x$specMod)
	invisible(x)
    }

print.summary.gel <- function(x, digits = 5, ...)
    {
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	if (is.null(x$CGEL))
            cat("Type of GEL: ", x$typeDesc, "\n")
	else
            cat("CGEL of type: ", x$typeDesc, " (alpha = ", x$CGEL, ")\n")
        
	if (!is.null(x$smooth))
            {
		cat("Kernel: ", x$smooth$kernel," (bw=", x$smooth$bw,")\n\n")
            }else {
		cat("\n")
            }
        
	cat("Coefficients:\n")
	print.default(format(x$coefficients, digits = digits),
                      print.gap = 2, quote = FALSE)
        
        if (length(x$coefficients)<length(x$lambda))
            {
                cat("\nLambdas:\n")
                print.default(format(x$lambda, digits=digits),
                              print.gap = 2, quote = FALSE)
            } else {
                cat("\nNo table for Lambda when the model is just identified\n")
            }
        cat("\n", x$stest$ntest, "\n")
        print.default(format(x$stest$test, digits=digits),
                      print.gap = 2, quote = FALSE)
        cat("\n",x$specMod)
        cat("\nConvergence code for the coefficients: ", x$conv_par, "\n")
        cat("\nConvergence code for the lambdas: ", x$conv_lambda$convergence, "\n")
        
        invisible(x)
    }

summary.gel <- function(object, ...)
	{
	z <- object
	n <- nrow(z$gt)
	se_par <- sqrt(diag(z$vcov_par))
	par <- z$coefficients
	tval <- par/se_par

	se_parl <- sqrt(diag(z$vcov_lambda))
	lamb <- z$lambda
	tvall <- lamb/se_parl

	ans <- list(type = z$type, call = z$call)
	names(ans$type) <-"Type of GEL"
	
	ans$coefficients <- round(cbind(par, se_par, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)), 5)
	ans$lambda <- round(cbind(lamb,se_parl, tvall, 2 * pnorm(abs(tvall), lower.tail = FALSE)), 5)

    	dimnames(ans$coefficients) <- list(names(z$coefficients), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    	dimnames(ans$lambda) <- list(names(z$lambda), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$stest=specTest(z)

	if (z$type == "EL")
		ans$badrho <- z$badrho
	if (!is.null(z$weights))
		{
		ans$weights <- z$weights
		}
	ans$conv_par <- z$conv_par
	ans$conv_pt <- z$conv_pt
	ans$conv_moment <- cbind(z$conv_moment)
	ans$conv_lambda <- z$conv_lambda
	ans$CGEL <- z$CGEL
        ans$typeDesc <- z$typeDesc
        ans$specMod <- z$specMod
	if (!is.null(attr(object$dat,"smooth")))
		ans$smooth <- attr(object$dat,"smooth")
	names(ans$conv_pt) <- "Sum_of_pt"
	dimnames(ans$conv_moment) <- list(names(z$gt), "Sample_moment_with_pt")
	class(ans) <- "summary.gel"
	ans	
}

residuals.gel <- function(object, ...) 
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$residuals
	}

fitted.gel <- function(object, ...)
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$fitted.value
	}

formula.gel <- function(x, ...)
{
    if(is.null(x$terms))
	stop("The gel object was not created by a formula")
    else
	formula(x$terms)
}

estfun.gel <- function(x, ...)
  {
  stop("estfun is not yet available for gel objects")
  }

bread.gel <- function(x, ...)
  {
  stop("Bread is not yet available for gel objects")
  }


getImpProb <- function(object, ...)
    UseMethod("getImpProb")

getImpProb.gel <- function(object, posProb=TRUE, normalize=TRUE,
                           checkConv=FALSE, ...)
    {
        if (!normalize || (object$type == "CUE" && !posProb))
            {
                n <- NROW(object$gt)
                pt <- -.rho(object$gt, object$lambda, type = object$type,
                            derive = 1, k = object$k1/object$k2)/n
                if (object$type == "CUE" && posProb)
                    {
                        eps <- -length(pt)*min(min(pt),0)
                        pt <- (pt+eps/length(pt))/(1+eps)
                    }
                if (normalize)
                    pt <- pt/sum(pt)
            } else {
                pt <- object$pt
            }
        if (checkConv)
            attr(pt, "convergence") <- list(pt=sum(pt),
                                            ptgt=colSums(pt*as.matrix(object$gt)))
        pt
    }


        
