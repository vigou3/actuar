### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Main interface to credibility model fitting functions.
###
### AUTHORS: Louis-Philippe Pouliot, Tommy Ouellet,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>.

cm <- function(formula, data, ratios, weights, subset,
               regformula = NULL, regdata, adj.intercept = FALSE,
               method = c("Buhlmann-Gisler", "Ohlsson", "iterative"),
               likelihood, ...,
               tol = sqrt(.Machine$double.eps), maxit = 100,
               echo = FALSE)
{
    Call <- match.call(expand.dots = TRUE)

    ## Catch the pure bayesian special case.
    if (formula == "bayes")
    {
        if (missing(data) || length(data) == 0L)
            data <- NULL
        res <- bayes(data, likelihood, ...)
        class(res) <- c("cm", class(res))
        attr(res, "call") <- Call
        return(res)
    }

    ## === MODEL ANALYSIS ===
    ##
    ## Decompose the formula giving the portfolio structure. Attribute
    ## "order" gives the interaction level of each term in the
    ## formula. In hierarchical structures, each term should represent
    ## a different level, hence there should not be any duplicates in
    ## this attribute. The column names in 'data' containing the
    ## portfolio structure can be obtained from the rownames of the
    ## matrix in attribute "factors".
    ##
    ## Note that the very last level, the data, is not taken into
    ## account here.
    tf <- terms(formula)
    level.numbers <- attr(tf, "order")           # level IDs
    level.names <- rownames(attr(tf, "factors")) # level names
    nlevels <- length(level.names)               # number of levels

    ## Sanity checks
    ##
    ## 1. only hierarchical interactions are allowed in 'formula';
    ## 2. hierarchical regression models are not supported.
    ## 3. if 'ratios' is missing, all columns of 'data' are taken to
    ##    be ratios, so 'weights' should also be missing;
    ##
    if (any(duplicated(level.numbers)))
        stop("unsupported interactions in 'formula'")
    if (nlevels > 1 && !is.null(regformula))
        stop("hierarchical regression models not supported")
    if (missing(ratios) & !missing(weights))
        stop("ratios have to be supplied if weights are")

    ## === DATA EXTRACTION ===
    ##
    ## 'data' is split into three matrices: one for the portfolio
    ## structure, one for the ratios and one for the weights. They are
    ## obtained via calls to subset() built from this function's
    ## call. That way, arguments 'ratios', 'weights' and 'subset' are
    ## not evaluated before being passed to subset(). Argument
    ## matching is as follows:
    ##
    ##   Argument of cm()     Argument of subset()
    ##   ================     ====================
    ##     data                 x
    ##     ratios               select
    ##     weights              select
    ##     subset               subset
    ##
    ## Positions of the arguments that will be needed.
    m <- match(c("data", "ratios", "weights", "subset"), names(Call), 0)

    ## Extraction of the portfolio structure. Arguments 'data' and
    ## 'subset' are passed to subset().
    cl <- Call[c(1, m[c(1, 4)])]        # use data and subset only
    cl[[1]] <- as.name("subset")        # change function name
    names(cl)[2] <- "x"                 # argument matching
    cl$select <- level.names            # add argument 'select'
    levs <- eval(cl, parent.frame())    # extraction

    ## Object 'levs' is a data frame or matrix with as many colums as
    ## there are levels in the model (still notwithstanding the data
    ## level). Rows contain nodes identifiers which can be
    ## anything. For calculations, these identifiers are converted
    ## into simple subscripts (i, j, k, ...) as used in mathematical
    ## notation.
    ##
    ## Note that 'apply' will coerce to a matrix.
    ilevs <- apply(levs, 2, function(x) as.integer(factor(x)))

    ## Extraction of the ratios. If argument 'ratios' is missing, then
    ## use all columns of 'data' except those of the portfolio
    ## structure.
    cl$select <-
        if (missing(ratios))
            setdiff(colnames(data), level.names)
        else
            Call[[m[2]]]
    ratios <- as.matrix(eval(cl, parent.frame())) # ratios as matrix

    ## Creation of a weight matrix. Extract from data if argument
    ## 'weights' is specified, otherwise create a matrix of ones. For
    ## extraction, the only change from ratio extraction is the
    ## content of element "select" of the call.
    weights <-
        if (missing(weights))
        {
            if (any(is.na(ratios)))
                stop("missing ratios not allowed when weights are not supplied")
            array(1, dim(ratios))       # matrix of ones
        }
        else
        {
            cl$select <- Call[[m[3]]]
            as.matrix(eval(cl, parent.frame())) # weights as matrix
        }

    ## == DISPATCH TO APPROPRIATE CALCULATION FUNCTION ==
    ##
    ## Buhlmann-Straub models are handled by bstraub(), regression
    ## models by hache() and hierarchical models by hierarc().
    if (nlevels < 2)                    # one-dimensional model
    {
        ## One-dimensional models accept only "unbiased" and
        ## "iterative" for argument 'method'.
        method <- match.arg(method)
        if (method == "Buhlmann-Gisler" || method == "Ohlsson")
            method <- "unbiased"

        if (is.null(regformula))        # Buhlmann-Straub
        {
            res <- bstraub(ratios, weights, method = method,
                           tol = tol, maxit = maxit, echo = echo)
        }
        else                            # Hachemeister
        {
            ## If regression model is actually empty or has only an
            ## intercept, call bstraub().
            trf <- terms(regformula)
            res <-
                if (length(attr(trf, "factors")) == 0)
                {
                    warning("empty regression model; fitting with Buhlmann-Straub's model")
                    bstraub(ratios, weights, method = method,
                            tol = tol, maxit = maxit, echo = echo)
                }
                else
                    hache(ratios, weights, regformula, regdata,
                          adj.intercept = adj.intercept,
                          method = method, tol = tol,
                          maxit = maxit, echo = echo)
        }

        ## Add missing quantities to results.
        res$classification <- levs
        res$ordering <- list(seq_along(levs))
    }
    else                                # hierarchical model
    {
        ## Computations with auxiliary function.
        res <- hierarc(ratios, weights, classification = ilevs,
                       method = method, tol = tol, maxit = maxit,
                       echo = echo)

        ## Put back original level names into the object
        res$classification <- levs
    }

    ## Transfer level names to lists
    names(res$means) <- names(res$weights) <- c("portfolio", level.names)
    names(res$unbiased) <- if (!is.null(res$unbiased)) names(res$means)
    names(res$iterative) <- if (!is.null(res$iterative)) names(res$means)
    names(res$nodes) <- names(res$ordering) <- level.names
    if (is.list(res$cred))
        names(res$cred) <- level.names

    ## Results
    class(res) <- c("cm", class(res))
    attr(res, "call") <- Call
    res
}

predict.cm <- function(object, levels = NULL, newdata, ...)
{
    ## Convert the character 'levels' argument into numeric and pass
    ## to next method.
    level.names <- names(object$nodes)

    levels <-
        if (is.null(levels))
            seq_along(level.names)
        else
            pmatch(levels, level.names)
    if (any(is.na(levels)))
        stop("invalid level name")
    NextMethod()
}

print.cm <- function(x, ...)
{
    chkDots(...)                        # method does not use '...'
    nlevels <- length(x$nodes)
    level.names <- names(x$nodes)
    b <- if (is.null(x$iterative)) x$unbiased else x$iterative

    cat("Call:\n",
        paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")

    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:", x$means[[1]], "\n")
    for (i in seq.int(nlevels))
    {
        if (i == 1L)
        {
            ## Treat the Hachemeister model separately since in this
            ## case the variance components vector is a list, with the
            ## first element a matrix. (Note that since a matrix with
            ## empty column names is printed to the screen, there will
            ## be a blank line in the display. Hence the inserted
            ## newline in the 'else' case.)
            if (attr(x, "model") == "regression")
            {
                m <- b[[1]]
                dimnames(m) <- list(c(paste("  Between", level.names[i], "variance: "),
                                      rep("", nrow(m) - 1)),
                                    rep("", ncol(m)))
                print(m)
            }
            else
                cat("\n  Between", level.names[i], "variance:",
                    b[i], "\n")
        }
        else
            cat("  Within ", level.names[i - 1],
                "/Between ", level.names[i], " variance: ",
                b[i], "\n", sep = "")
    }
    cat("  Within", level.names[nlevels], "variance:",
        b[[nlevels + 1]], "\n", fill = TRUE)
    invisible(x)
}

summary.cm <- function(object, levels = NULL, newdata, ...)
{
    nlevels <- length(object$nodes)

    if (nlevels == 1L)
    {
        ## Single level cases (Buhlmann-Straub and Hachemeister):
        ## return the object with the following modifications: put
        ## credibility factors into a list and add a list of the
        ## credibility premiums.
        object$premiums <- list(predict(object, newdata = newdata))
        object$cred <- list(object$cred)
        class(object) = c("summary.cm", class(object))
    }
    else
    {
        ## Multi-level case (hierarchical): select result of the
        ## appropriate level(s).
        plevs <-
            if (is.null(levels))
                seq_along(names(object$nodes))
            else
                pmatch(levels, names(object$nodes))
        if (any(is.na(plevs)))
            stop("invalid level name")

        object$premiums <- predict(object, levels) # new element
        object$means <- object$means[c(1, plevs + 1)]
        object$weights <- object$weights[c(1, plevs + 1)]
        object$unbiased <- object$unbiased[sort(unique(c(plevs, plevs + 1)))]
        object$iterative <- object$iterative[sort(unique(c(plevs, plevs + 1)))]
        object$cred <- object$cred[plevs]
        object$classification <- object$classification[, seq.int(max(plevs)), drop = FALSE]
        object$nodes <- object$nodes[plevs]
        class(object) <- c("summary.cm", class(object))
    }
    structure(object, ...)     # attach additional attributes in '...'
}

print.summary.cm <- function(x, ...)
{
    nlevels <- length(x$nodes)
    level.names <- names(x$nodes)
    NextMethod()                        # print.cm()
    cat("Detailed premiums\n\n")
    for (i in seq.int(nlevels))
    {
        ## Print a "section title" only if there is more than one
        ## level. (Provision introduced in v2.3.0; before the title
        ## was always printed.)
        if (nlevels > 1L)
            cat("  Level:", level.names[i], "\n")

        ## There are no level names in the linear Bayes case, so we
        ## skip this column in the results.
        if (is.null(level.names))
            levs <- NULL
        else
        {
            level.id <- match(level.names[i], colnames(x$classification))
            levs <- x$classification[, seq.int(level.id), drop = FALSE]
            m <- duplicated(levs)
        }

        if (attr(x, "model") == "regression")
        {
            ## Hachemeister model: results contain matrices
            y <- cbind(" ",
                       as.vector(format(x$means[[i + 1L]], ...)),
                       as.vector(apply(format(x$cred[[i]], ...), c(1L, 3L),
                                       paste, collapse = " ")),
                       as.vector(format(sapply(x$adj.models, coef), ...)),
                       " ")
            y[seq(1, nrow(y), dim(x$cred[[i]])[1]), c(1L, 5L)] <-
                c(levs[!m, , drop = FALSE], format(x$premiums[[i]], ...))
            colnames(y) <- c(colnames(levs),
                             "Indiv. coef.", "Cred. matrix",
                             "Adj. coef.", "Cred. premium")
        }
        else if (is.null(levs))
        {
            ## Linear Bayes model: simplified results with no level
            ## column
            y <- cbind(format(x$means[[i + 1L]], ...),
                       format(x$weights[[i + 1L]], ...),
                       format(x$cred[[i]], ...),
                       format(x$premiums[[i]], ...))
            colnames(y) <- c("Indiv. mean", "Weight",
                             "Cred. factor", "Bayes premium")
        }
        else
        {
            ## All other models
            y <- cbind(as.matrix(levs[!m, , drop = FALSE]),
                       format(x$means[[i + 1L]], ...),
                       format(x$weights[[i + 1L]], ...),
                       format(x$cred[[i]], ...),
                       format(x$premiums[[i]], ...))
            colnames(y) <- c(colnames(levs),
                             "Indiv. mean", "Weight",
                             "Cred. factor", "Cred. premium")
        }
        rownames(y) <- rep("   ", nrow(y))
        print(y, quote = FALSE, right = FALSE, ...)
        cat("\n")
    }
    invisible(x)
}
