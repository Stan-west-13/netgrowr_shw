#' Compare two models with likelihood ratio test
#'
#' @param full list object with class \code{\link[netgrowr]{mle_network_growth}} associated with the full model.
#' @param restricted list object with class \code{\link[netgrowr]{mle_network_growth}} associated with the restricted (i.e., nested) model.
#' @return A named vector with the following values:
#' \item{p}{One-tailed probability that the likelihood ratio is sampled from a
#' chi-square distribution with \code{df} degrees of freedom centered on zero or
#' less.}
#' \item{df}{The degrees of freedom of the statistical test; the difference in
#' the number of parameters between the nested models.}
#' \item{nobs}{The number of examples used to fit the models; the number of fitted.values for each model (which must be the same).}
#' \item{L0}{The negative log-likelihood of the restricted model.}
#' \item{L1}{The negative log-likelihood of the full model.}
#' \item{chisq}{The chi-square value of the likelihood ratio.}
#' \item{BIC}{The Bayesian Information Criterion for the model comparison, based on the chi-squared likelihood ratio.}
#'
#' @export
model_comparison <- function(full, restricted) {
    assertthat::are_equal(full$nobs, restricted$nobs)
    df <- length(coef(full)) - length(coef(restricted))
    chisq <- likelihood_ratio(full$negLogLik, restricted$negLogLik)
    pval <- pchisq(chisq, df = df, lower.tail = FALSE)
    BIC <- bayesian_information_criterion(chisq, df, full$nobs)
    return(c(
        df = df,
        nobs = full$nobs,
        L0 = restricted$negLogLik,
        L1 = full$negLogLik,
        chisq = chisq,
        p = pval,
        BIC = BIC
    ))
}

#' likelihood ratio chi-square value
#'
#' Compares the negative log-likelihood between a pair of nested models,
#' producing a value that follows the chi-square distribution.
#'
#' @param neg_loglik_full Negative log-likelihood of the full model
#' @param neg_loglik_restricted Negative log-likelihood of the restricted model
#' @return Chi-squared value
#'
likelihood_ratio <- function(neg_loglik_full, neg_loglik_restricted) {
    return(2 * (neg_loglik_restricted - neg_loglik_full))
}

#' Bayesian information criterion (BIC)
#'
#' @param chisq The chi-square value of the likelihood ratio.
#' @param df The degrees of freedom for the likelihood ratio test.
#' @param n The number of observations the models were fit to.
#' @return A real scalar
#'
bayesian_information_criterion <- function(chisq, df, nobs) {
    return(chisq - (df * log(nobs)))
}
