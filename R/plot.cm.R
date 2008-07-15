### ===== actuar: an R package for Actuarial Science =====
### File trunk/R/plot.cm.R
###
### Plot an object of class "cm" and display the homogeneity of
### the portfolio showing the ratio between the within variance and
### the between variance. Strongly based on file plot.lm.R
###
### AUTHORS: Xavier Milhaud, Vincent Goulet
### <vincent.goulet@act.ulaval.ca>

plot.cm <-
function (x, which = c(1,2,3),
	  caption = c("Ratio within/between variance in Buhlmann-Straub model",
          "Ratio within/between variance in Hachemeister model",
          "Ratio within/between variance ina hierarchical model"),
          contractNo = NULL,
          method = c("iterative", "unbiased"),
	  panel = if(add.smooth) panel.smooth else points,
          main = "", ...,
          label.pos = c(4,2), cex.caption = 1)
{
    if (!inherits(x, "cm"))
	stop("use only with \"cm\" objects")
    if (!is.numeric(which) || (any(which) < 1) ||(any(which) > 3))
	stop("'which' must be one of the three models required")

    require(graphics)
    show <- rep(FALSE, 3)
    show[which] <- TRUE

    ##---------- Do the individual plots : ----------
    if (show[1] || show[3])   ## Buhlmann-Straub or hierarc models
    {
        method <- match.arg(method)
        if (method == "iterative")
            mat <- rbind(x$iterative[[1]][1,1], x$iterative[[2]])
        else
            mat <- rbind(x$unbiased[[1]][1,1], x$unbiased[[2]])
        ## draw a barplot with extreme values equals to within and
        ## between variance. The larger it is, the less homogeneous
        ## the portfolio is.
        boxplot(as.data.frame(mat), col = "grey", main = caption[2],
                xlab = "portfolio", ylab = "variance")
    }
    if (show[2])    ## Hachemeister model
    {
        ## Credibility, individual and collective regression lines
        mu_ind <- mu_coll <- mu_cred <- numeric(ncol(weights))
        for (i in 1:ncol(weights))
        {
            mu_ind[i] <- x$means[[2]][1, contractNo] + i * x$means[[2]][2, contractNo]
            mu_coll[i] <- x$means[[1]][1] + i * x$means[[1]][2]
            mu_cred[i] <- x$adj.models[[contractNo]]$coefficients[1] +
                i * x$adj.models[[contractNo]]$coefficients[1]
        }
        ylim = range(c(mu_ind, mu_coll))
        plot(1:ncol(weights), ratios[contractNo,], type = "p",
             ylim = extendrange(r = ylim, f = 0.05),
             xlab = "contract", ylab = "premiums")
        lines(1:ncol(weights), mu_ind, col = "red")
        lines(1:ncol(weights), mu_coll, col = "blue")
        lines(1:ncol(weights), mu_cred, col = "green")
        mtext(caption[2], 3, 0.25, cex = cex.caption)
        ## or SEE boxplot(formula,...) ??
    }
}
