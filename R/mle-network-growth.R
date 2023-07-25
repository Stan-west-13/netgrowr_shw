#' Fit model of network growth using Maximum Likelihood Estimation (MLE)
#'
#' Objective is to maximize the predicted probability of nodes being added to
#' the network at the point in time when they are known to have been acquired.
#'
#' @param formula A formula; see details
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @param beta A vector of weights (i.e, the model), or a matrix where each
#'   column is a different model.
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details.
#' @param label_with A string referencing a variable in \code{data} that
#'   provides node-labels. Labels may repeat across time steps specified by
#'   \code{split_by}.
#' @param B0 An initial estimate of the model parameters. If \code{NULL}, an
#'   initial model will be estimated through a brute force optimization.
#'   [default: NULL]
#' @return  A list as of class \code{mle_network_growth} with the following fields:
#' \item{coefficients}{Best fit model parameters}
#' \item{negLogLik}{The model negative log-likelihood.}
#' \item{fitted.values}{A vector of estimated probabilities for each node acquired over all timepoints.}
#' \item{nobs}{The number of fitted values}
#' \item{call}{The function call that generated the model}
#' \item{model}{The model frame referenced when fitting the model.}
#' \item{deviance}{The likelihood ratio between a fully saturated model and the current model}
#' \item{counts}{A two-element integer vector giving the number of calls to objective function and gradient function respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to the objective function to compute a finite-difference approximation to the gradient. This is returned by \code{\link[stats]{optim}}.}
#' \item{convergence}{An integer code:
#'     \describe{
#'         \item{0}{successful completion (which is always the case for "SANN" and "Brent").}
#'         \item{1}{indicates that the iteration limit maxit had been reached.}
#'         \item{10}{indicates degeneracy of the Nelder–Mead simplex.}
#'         \item{51}{indicates a warning from the "L-BFGS-B" method; see component message for further details.}
#'         \item{52}{indicates an error from the "L-BFGS-B" method; see component message for further details.}
#'     }
#' }
#' \item{message}{A character string giving any additional information returned by the optimizer, or NULL.}
#' \item{hessian}{Only included if argument \code{hessian} is true. A symmetric matrix giving an estimate of the Hessian at the solution found. Note that this is the Hessian of the unconstrained problem even if the box constraints are active.}
#'
#' @details
#' The formula must reference a variable that encodes the time-point at which
#' each node is acquired as the dependent variable. If \code{data}
#' contains the data for multiple time steps, then this value should probably be
#' the same for all points in time at which the node is referenced (although the
#' function does not enforce this). The independent variables specify the set of
#' variables that will be used to predict, at each point in time, whether the
#' node is acquired or not.
#'
#' The \code{split_by} variable is used to split \code{data} by time-step, so
#' that probabilities can be estimated independently at discrete points in time.
#' The assumption is that a node will be acquired once, at a single point in
#' time.
#'
#' The function will return a vector if a beta was a vector or one-column
#' matrix. Otherwise, the return value will be a matrix, with a row for each
#' node and a column for each model specified by \code{beta}. Regardless, the
#' returned values are probabilities for nodes that are acquired, as predicted
#' at the point in time the node was acquired. Probabilities are computed as a
#' ratio of strengths.
#'
#' @seealso \code{\link[netgrowr]{probability_node_added}}, \code{\link[netgrowr]{ratio_of_strengths}}
#'
#' @export
mle_network_growth <- function(formula, data, split_by, label_with, B0 = NULL) {
    if (is.null(B0)) {
        B0 <- bruteforce_optim(formula, data, split_by = split_by, R = 1e4)
    }
    M <- gradient_optim(formula, data, B0, split_by = split_by)
    names(M)[names(M) == "value"] <- "negLogLik"
    names(M)[names(M) == "par"] <- "coefficients"
    names(M$coefficients) <- if (attr(terms(formula), "intercept")) {
        c("(intercept)", attr(terms(formula), "term.labels"))
    } else {
         attr(terms(formula), "term.labels")
    }
    p <- probability_node_added(coef(M), formula, data, split_by = split_by, label_with = label_with)
    M$fitted.values <- p
    M$nobs <- length(p)
    M$call <- match.call()
    M$model <- model_frame_with_split_and_labels(formula, split_by, label_with, data)
    learned <- M$model[[1]] == M$model[[split_by]]
    p_star <- unlist(lapply(
        split(learned[learned], M$model[[split_by]][learned]),
        function(x) x / sum(x)
    ))
    M$deviance <- likelihood_ratio(-loglikelihood(p_star), -loglikelihood(p))
    attr(M, "class") <- "mle_network_growth"
    return(M)
}

#' @method print mle_network_growth
#' @export
print.mle_network_growth <- function(x, digits = max(3L, getOption("digits") - 3L)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                           collapse = "\n"), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L,
                      quote = FALSE)
    } else {
        cat("No coefficients\n")
    }
    cat("\n")
    performance <- c(
        "negative log-likelihood" = x$negLogLik,
        "deviance" = x$deviance
    )
    cat("Performance:\n")
    print.default(format(performance, digits = digits), print.gap = 2L,
                  quote = FALSE)
    cat("\n")
    cat("Number of observations:", x$nobs, "\n")
    cat("\n")
    invisible(x)
}

#' Log likelihood
#'
#' @param p A vector or matrix of probabilities; see details.
#' @return If \code{p} is a matrix, a vector of log likelihoods as long as \code{ncol(p)}, otherwise a single value.
#'
#' @details
#' If \code{p} is a matrix, each column represents the probabilities associated with a different model.
#'
loglikelihood <- function(p) {
    if (is.matrix(p)) {
        return(colSums(log(p), na.rm = TRUE))
    } else {
        return(sum(log(p), na.rm = TRUE))
    }
}

#' Brute force optimization
#'
#' Evaluate a large number of randomly generated models. Useful for obtaining an
#' initial model to begin a gradient optimization.
#'
#' @param formula A formula; see details
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details.
#' @param R The number of random models to evaluate.
#' @return A vector expressing the best model parameters discovered in the brute
#'   force search.
#'
#' @details
#' The formula must reference a variable that encodes the time-point at which
#' each node is acquired as the dependent variable. If \code{data}
#' contains the data for multiple time steps, then this value should probably be
#' the same for all points in time at which the node is referenced (although the
#' function does not enforce this). The independent variables specify the set of
#' variables that will be used to predict, at each point in time, whether the
#' node is acquired or not.
#'
#' The \code{split_by} variable is used to split \code{data} by time-step, so
#' that probabilities can be estimated independently at discrete points in time.
#' The assumption is that a node will be acquired once, at a single point in
#' time.
bruteforce_optim <- function(formula, data, split_by, R = 1e4) {
    nbeta <- nterms(formula)
    Brandom <- matrix(runif(nbeta * R, min = -1.0, max = 1.0), nrow = nbeta, ncol = R)
    p <- probability_node_added(Brandom, formula = formula, data = data, split_by = split_by)
    return(Brandom[ , which.max(loglikelihood(p))])
}

#' Gradient optimization of the negative log likelihood
#'
#' Nelder-Mead optimization using the \code{\link[stats]{optim}} function.
#'
#' @param formula A formula; see details
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @param B0 Initial values for the model parameters.
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details.
#' @param control See documentation for \code{\link[stats]{optim}}.
#' @return See documentation for \code{\link[stats]{optim}}.
#'
#' @details
#' The optimization minimizes the negative log likelihood associated with the
#' estimated probabilities for each nodes at the time it is added to the
#' network.
#'
#' The formula must reference a variable that encodes the time-point at which
#' each node is acquired as the dependent variable. If \code{data}
#' contains the data for multiple time steps, then this value should probably be
#' the same for all points in time at which the node is referenced (although the
#' function does not enforce this). The independent variables specify the set of
#' variables that will be used to predict, at each point in time, whether the
#' node is acquired or not.
#'
#' The \code{split_by} variable is used to split \code{data} by time-step, so
#' that probabilities can be estimated independently at discrete points in time.
#' The assumption is that a node will be acquired once, at a single point in
#' time.
gradient_optim <- function(formula, data, B0, split_by, control = list(maxit = 1e3)) {
    nLL <- function(beta, formula, data) {
        p <- probability_node_added(beta = beta, formula = formula, data = data, split_by = split_by)
        return(-loglikelihood(p))
    }
    method <- if (length(B0) == 1) 'Brent' else 'Nelder-Mead'
    lower <- if (length(B0) == 1) -1 else -Inf
    upper <- if (length(B0) == 1) 1 else Inf
    return(optim(par = B0, fn = nLL, method = method, control = control,
                 formula = formula, lower = lower, upper = upper, data = data))
}

#' Compute ratio of strengths
#'
#' Standardize the evidence for the targets by dividing by the sum of the
#' evidence across targets and distractors.
#'
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @param formula A formula; see details
#' @param beta A vector of weights (i.e, the model), or a matrix where each
#'   column is a different model.
#' @return A matrix of probabilities associated with nodes that were
#'   acquired. There will be a column for each model specified by \code{beta}.
#'
#' @details
#' The formula must reference a variable that encodes the time-point at which
#' each node is acquired as the dependent variable. If \code{data}
#' contains the data for multiple time steps, then this value should probably be
#' the same for all points in time at which the node is referenced (although the
#' function does not enforce this). The independent variables specify the set of
#' variables that will be used to predict, at each point in time, whether the
#' node is acquired or not.
#'
ratio_of_strengths <- function(data, formula, beta) {
    learned <- data$learned
    unknown <- data$unknown
    d <- model.matrix(formula, data, na.action = "na.fail")
    x <- exp(d %*% as.matrix(beta))
    # x is a matrix with a row for each item and a column for each set of betas.
    # colSums will produce a vector with an element for each set of betas.
    # Transposing x facilitates efficient ratio computation.
    # The resulting matrix is transposed to get back to the original orientation.
    return(t(t(x[learned, , drop = FALSE]) / colSums(x[unknown, , drop = FALSE])))
}

#' Generate predicted acquisition probabilities for network nodes
#'
#' Compute the ratio of strengths to evaluate the model evidence for acquiring a
#' node at the time it was acquired. The \code{split_by} parameter
#' specifies a variable in \code{data} that specifies different points in time.
#'
#' @param beta A vector of weights (i.e, the model), or a matrix where each
#'   column is a different model.
#' @param formula A formula; see details
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details.
#' @param label_with A string referencing a variable in \code{data} that
#'   provides node-labels. Labels may repeat across time steps specified by
#'   \code{split_by}.
#' @return A numeric vector or matrix; see details.
#'
#' @details
#' The formula must reference a variable that encodes the time-point at which
#' each node is acquired as the dependent variable. If \code{data}
#' contains the data for multiple time steps, then this value should probably be
#' the same for all points in time at which the node is referenced (although the
#' function does not enforce this). The independent variables specify the set of
#' variables that will be used to predict, at each point in time, whether the
#' node is acquired or not.
#'
#' The \code{split_by} variable is used to split \code{data} by time-step, so
#' that probabilities can be estimated independently at discrete points in time.
#' The assumption is that a node will be acquired once, at a single point in
#' time.
#'
#' The function will return a vector if a beta was a vector or one-column
#' matrix. Otherwise, the return value will be a matrix, with a row for each
#' node and a column for each model specified by \code{beta}. Regardless, the
#' returned values are probabilities for nodes that are acquired, as predicted
#' at the point in time the node was acquired. Probabilities are computed as a
#' ratio of strengths.
#'
#' @seealso \code{\link[netgrowr]{ratio_of_strengths}}
#'
probability_node_added <- function(beta, formula, data, split_by, label_with = NULL) {
    d <- na.omit(get_all_vars_with_split_and_labels(formula, split_by, label_with, data))
    d$learned <- d[[1]] == d[[split_by]]
    d$unknown <- d[[1]] >= d[[split_by]]
    X <- split(d, d[[split_by]])
    if (!is.null(label_with)) {
        X <- lapply(X, function(x) {
                rownames(x) <- x[[label_with]]
                return(x)
            })
    }
    p <- do.call("rbind", lapply(X, netgrowr:::ratio_of_strengths, beta = beta, formula = formula))
    if (ncol(p) == 1) {
        labs <- rownames(p)
        p <- as.vector(p)
        names(p) <- labs
    }
    return(p)
}

#' Count terms (independent variables) implied by a function
#'
#' Unless the formula explicitly states that an intercept term is not desired,
#' the returned count includes an intercept term.
#'
#' @param formula a formula
#' @return The number of terms the formula implies.
nterms <- function(formula) {
    x <- terms(formula)
    return(length(attr(x, "term.labels")) + attr(x, "intercept"))
}

#' Return model frame with split and and label variables
#'
#' Written to handle the case where either split or labels are not provided.
#'
#' @param formula a formula
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details for \code{\link[netgrowr]{mle_network_growth}}.
#' @param label_with A string referencing a variable in \code{data} that
#'   provides node-labels. Labels may repeat across time steps specified by
#'   \code{split_by}.
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @return The model frame
model_frame_with_split_and_labels <- function(formula, split_by, label_with, data) {
    f <- formula_with_split_and_labels(formula, split_by, label_with)
    return(model.frame(f, data))
}

#' Return variables from fromula with split and and label variables
#'
#' Written to handle the case where either split or labels are not provided.
#'
#' @param formula a formula
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details for \code{\link[netgrowr]{mle_network_growth}}.
#' @param label_with A string referencing a variable in \code{data} that
#'   provides node-labels. Labels may repeat across time steps specified by
#'   \code{split_by}.
#' @param data A data frame containing, at least, the variables implied by the
#'   formula.
#' @return The model frame, without transformations implied in the formula.
get_all_vars_with_split_and_labels <- function(formula, split_by, label_with, data) {
    f <- formula_with_split_and_labels(formula, split_by, label_with)
    return(get_all_vars(f, data))
}

#' Return model frame with split and and label variables
#'
#' Written to handle the case where either split or labels are not provided.
#'
#' @param formula a formula
#' @param split_by A string referencing a variable in \code{data} that encodes
#'   time steps; see details for \code{\link[netgrowr]{mle_network_growth}}.
#' @param label_with A string referencing a variable in \code{data} that
#'   provides node-labels. Labels may repeat across time steps specified by
#'   \code{split_by}.
#' @return A formula with split and label variables appended.
formula_with_split_and_labels <- function(formula, split_by, label_with) {
    f <- formula
    if (!is.null(split_by))
        f <- update(f, paste("~.", split_by, sep = "+"))
    if (!is.null(label_with))
        f <- update(f, paste("~.", label_with, sep = "+"))
    return(f)
}
