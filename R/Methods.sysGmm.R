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

summary.sysGmm <- function(object, ...)
    {
	z <- object
	se <- sqrt(diag(z$vcov))
        k <- attr(z$dat, "k")
        if (attr(z$dat, "sysInfo")$commonCoef)
            {
                seList <- rep(list(se), length(z$dat))
            } else {
                seList <- list()
                for (i in 1:length(z$dat))
                    {
                        seList[[i]] <- se[1:k[[i]]]
                        se <- se[-(1:k[[i]])]
                    }
            }
	par <- z$coefficients
	tval <- lapply(1:length(z$dat), function(i) par[[i]]/seList[[i]])
	ans <- list(met=z$met,kernel=z$kernel,algo=z$algo,call=z$call)
	names(ans$met) <- "GMM method"
	names(ans$kernel) <- "kernel for cov matrix"
        
	ans$coefficients <- lapply(1:length(z$dat), function(i) cbind(par[[i]],seList[[i]], tval[[i]], 2 * pnorm(abs(tval[[i]]), lower.tail = FALSE)))
        ans$stest <- specTest(z)
        ans$algoInfo <- z$algoInfo
        ans$initTheta <- object$initTheta

        for (i in 1:length(z$dat))
            {
            dimnames(ans$coefficients[[i]]) <- list(names(z$coefficients[[i]]), 
                                                    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
            names(ans$initTheta[[i]]) <- names(z$coefficients[[i]])
        }
        ans$specMod <- object$specMod
	ans$bw <- attr(object$w0,"Spec")$bw
	ans$weights <- attr(object$w0,"Spec")$weights
        ans$Sysnames <- names(z$dat)
        ans$met <- attr(object$dat, "sysInfo")$typeDesc
	if(object$infVcov != "HAC")
            ans$kernel <- NULL
	class(ans) <- "summary.sysGmm"
	ans
    }

print.summary.sysGmm <- function(x, digits = 5, ...)
    {
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
        cat("Method\n", x$met,"\n\n")
        cat("\n")
	if( !is.null(x$kernel))
            {
		cat("Kernel: ", x$kernel)
		if (!is.null(x$bw))
                    cat("(with bw = ", round(x$bw,5),")\n\n")
		else
                    cat("\n\n")	
            }
	cat("Coefficients:\n")
        m <- length(x$coefficients)
        for (i in 1:m)
            {
                cat(x$Sysnames[i], "\n")
                cat("#########\n")
                #print.default(format(x$coefficients[[i]], digits=digits),
                #              print.gap = 2, quote = FALSE)
                printCoefmat(x$coefficients[[i]], digits=digits, ...)
                cat("\n")
            }
	cat(x$stest$ntest,"\n")
	print.default(format(x$stest$test, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	if(!is.null(x$initTheta))
            {
		cat("Initial values of the coefficients\n")
                for (i in 1:m)
                    {
                        cat(x$Sysnames[i], "\n")
                        print(x$initTheta[[i]])
                    }
		cat("\n")
            }
        cat(x$specMod)
	if(!is.null(x$algoInfo))
            {	
		cat("#############\n")
 		cat("Information related to the numerical optimization\n")
            }
	if(!is.null(x$algoInfo$convergence))
            cat("Convergence code = ", x$algoInfo$convergence,"\n")
	if(!is.null(x$algoInfo$counts))
            {	
		cat("Function eval. = ",x$algoInfo$counts[1],"\n")
		cat("Gradian eval. = ",x$algoInfo$counts[2],"\n")
            }	
	if(!is.null(x$algoInfo$message))
            cat("Message: ",x$algoInfo$message,"\n")
	invisible(x)
    }

print.sysGmm <- function(x, digits=5, ...)
    {
	cat("Method\n", attr(x$dat, "sysInfo")$typeDesc,"\n\n")
	cat("Objective function value: ",x$objective,"\n\n")
        for (i in 1:length(x$coefficients))
            {
                cat(names(x$dat)[[i]], ": \n")
                print.default(format(coef(x)[[i]], digits=digits),
                              print.gap = 2, quote = FALSE)
            }
	cat("\n")
	if(!is.null(x$algoInfo$convergence))
            cat("Convergence code = ", x$algoInfo$convergence,"\n")
	cat(x$specMod)
	invisible(x)
    }





		


