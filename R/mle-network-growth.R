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
#'   provides node-labels.
#' @param B0 An initial estimate of the model parameters. If \code{NULL}, an
#'   initial model will be estimated through a brute force optimization.
#'   [default: NULL]
#' @return  A list as returned by \code{\link[stats]{optim}}, with the following additional fields:
#' \item{p}{A vector of estimated probabilities for each node acquired over all timepoints.}
#' \item{nobs}{The number of nodes for which a probability wsa estimated.}
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
    p <- probability_node_added(M$par, formula, data, split_by = split_by, label_with = label_with)
    M$p <- p
    M$nobs <- length(p)
    return(M)
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
    d <- model.matrix(formula, data, na.action = "na.omit")
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
#'   provides node-labels.
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
    f <- update(formula, trimws(paste("~.", split_by, label_with, sep = "+"), which = "right", whitespace = "\\+"))
    d <- na.omit(get_all_vars(f, data))
    d$learned <- d[[1]] == d[[split_by]]
    d$unknown <- d[[1]] >= d[[split_by]]
    X <- split(d, d[[split_by]])
    if (!is.null(label_with)) {
        X <- lapply(X, function(x) {
                rownames(x) <- x[[label_with]]
                return(x)
            })
    }
    f <- update
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
