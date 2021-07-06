#' Compute growth values based on network structure
#'
#' @param net A binary adjacency matrix.
#' @param aoa_tbl A vector of known vertexes. Can be a logical vector.
#' @param weighted Logical. Should weights be taken into account?
#' @param alpha What proportion of the growth value should be influenced by the
#'   weight [1]? This is ignored if \code{weighted} is FALSE.
#' @return A vector of growth values for unknown words.
#'
#' Growth values are computed for the words represented by each column. If the
#' model is not expressed by a symmetric matrix, the column names are used to
#' filter the rows so that it is made to be symmetric when computing growth values
#' according to the lure of the associates and preferential attechment. These
#' models are concerned with only the words that can be learned---and in the
#' context of these growth values, only words that have AoA values can be learned.
#' So the columns correspond to words for which we have AoA values, and that
#' defines the scope of the growth model.
#'
#' The preferential acquisition growth model, unlike the other two, takes into
#' account the structure of the whole environment, not limited to only those words
#' that can be learned. In other words, the indegree of a word that can be learned
#' can be influenced by words that are in the environment but cannot be learned.
#' Thus, the model can be asymmetric, such that it has more rows than columns,
#' and this will change the pattern of growth predicted by preferential
#' acquisition.
#'
#' @seealso \code{\link[netgrowr]{PreferentialAttachment}}, \code{\link[netgrowr]{LureOfTheAssociates}}, \code{\link[netgrowr]{PreferentialAcquisition}}, \code{\link[netgrowr]{opsahl_weighted_degree}}
#'
#' @references
#' 1. Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32(3), 245-251. doi:https://doi.org/10.1016/j.socnet.2010.03.006
#'
#' @export
growth_values <- function(net, aoa_tbl, growth_models = "all", weighted = FALSE, alpha = 0.5) {
    aoa <- align_aoa_vec(net, aoa_tbl)
    words <- colnames(net)
    nwords <- length(words)
    months <- sort(unique(aoa[!is.na(aoa)]))[-1]
    nmonths <- length(months)
    aoa[is.na(aoa)] <- Inf

    known <- vapply(months, is_known, logical(nwords), aoa = aoa)
    learned <- vapply(months, is_learned, logical(nwords), aoa = aoa)

    if (length(growth_models) == 1 && growth_models == "all") {
        growth_models <- c("preferential_attachment", "lure_of_the_associates", "preferential_acquisition")
    }

    growth_values <- lapply(growth_models, function(g) {
        apply(known, 2, g, net = net, weighted = weighted, alpha = alpha)
    })
    names(growth_values) <- growth_models
    growth_valuesZ <- lapply(growth_values, scale)


    X <- lapply(growth_models, function(g) {
        data.frame(model = factor(g, growth_models),
                   word = factor(rep(words, nmonths), words),
                   month = rep(months, each = nwords),
                   aoa = rep(aoa, nmonths),
                   known = c(known),
                   learned = c(learned),
                   value = c(growth_values[[g]]),
                   zscore = c(growth_valuesZ[[g]]))
    })
    return(do.call("rbind", X))
}

#' Compute preferential attachment values
#'
#' @param net A binary adjacency matrix.
#' @param known A vector of known vertexes. Can be a logical vector.
#' @param weighted Logical. Should weights be taken into account?
#' @param alpha What proportion of the growth value should be influenced by the
#'   weight [1]? This is ignored if \code{weighted} is FALSE.
#' @return A vector of growth values for unknown words.
#'
#' Words that are flagged as \code{known} in the input will be associated with
#' \code{NA} in the output.
#'
#' Children acquire words that are related to high-degree words that they
#' already know.
#'
#' Growth value is the mean degree of the known words that the new word
#' attached to [2].
#'
#' @seealso \code{\link[netgrowr]{GrowthValue}}, \code{\link[netgrowr]{LureOfTheAssociates}}, \code{\link[netgrowr]{PreferentialAcquisition}}, \code{\link[netgrowr]{opsahl_weighted_degree}}, \code{\link[netgrowr]{opsahl_weighted_degree}}
#'
#' @references
#' 1. Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32(3), 245-251. doi:https://doi.org/10.1016/j.socnet.2010.03.006
#'
#' 2. Steyvers, M., & Tenenbaum, J. B. (2005). The large‐scale structure of semantic networks: Statistical analyses and a model of semantic growth. Cognitive science, 29(1), 41-78.
#'
#' @export
preferential_attachment <- function(net, known, weighted = FALSE, alpha = 0.5) {
    patt <- function(net, known) {
        X <- net[!known, known, drop = FALSE]
        nd <- as.matrix(colSums(net[known, known, drop = FALSE]))
        gv <- numeric(length(known))
        gv[!known] <- (X %*% nd) / rowSums(X > 0)
        gv[is.na(gv)] <- 0
        gv[is.infinite(gv)] <- 0
        gv[known] <- NA
        return(gv)
    }
    if (ncol(net) != nrow(net)) {
        net <- net[rownames(net) %in% colnames(net), , drop = FALSE]
    }
    diag(net) <- 0
    gv <- patt(net > 0, known)
    if (weighted) {
        gvw <- patt(net, known)
        gv <- opsahl_weighted_degree(gv, gvw, alpha)
    }
    names(gv) <- colnames(net)
    return(gv)
}

#' Compute lure of the associates
#'
#' @param net A binary adjacency matrix.
#' @param known A vector of known vertexes. Can be a logical vector.
#' @param weighted Logical. Should weights be taken into account?
#' @param alpha What proportion of the growth value should be influenced by the
#'   weight [1]? This is ignored if \code{weighted} is FALSE.
#' @return A vector of growth values for unknown words.
#'
#' Words that are flagged as \code{known} in the input will be associated with
#' \code{NA} in the output.
#'
#' Children acquire words that are associated with many words that are
#' already known.
#'
#' Growth value is the degree of the word with respect to links from known words at
#' the time of acquisition [2].
#'
#' @seealso \code{\link[netgrowr]{GrowthValue}}, \code{\link[netgrowr]{PreferentialAttachment}}, \code{\link[netgrowr]{PreferentialAcquisition}}, \code{\link[netgrowr]{opsahl_weighted_degree}}
#'
#' @references
#'   1. Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32(3), 245-251. doi:https://doi.org/10.1016/j.socnet.2010.03.006
#'
#'   2. Hills, T. T., Maouene, M., Maouene, J., Sheya, A., & Smith, L. (2009). Longitudinal analysis of early semantic networks: Preferential attachment or preferential acquisition? Psychological science, 20(6), 729-739. doi:10.1111/j.1467-9280.2009.02365.x
#'
#' @export
lure_of_the_associates <- function(net, known, weighted = FALSE, alpha = 0.5) {
    loa <- function(net, known) {
        X <- net[known, !known, drop = FALSE]
        gv <- numeric(length(known))
        gv[!known] <- colSums(X)
        gv[known] <- NA
        return(gv)
    }
    if (ncol(net) != nrow(net)) {
        net <- net[rownames(net) %in% colnames(net), , drop = FALSE]
    }
    diag(net) <- 0
    gv <- loa(net > 0, known)
    if (weighted) {
        gvw <- loa(net, known)
        gv <- opsahl_weighted_degree(gv, gvw, alpha)
    }
    names(gv) <- colnames(net)
    return(gv)
}

#' Compute preferential acquisition values
#'
#' @param net A binary adjacency matrix.
#' @param known A vector of known vertexes. Can be a logical vector.
#' @param weighted Logical. Should weights be taken into account?
#' @param alpha What proportion of the growth value should be influenced by the
#'   weight [1]? This is ignored if \code{weighted} is FALSE.
#' @return A vector of growth values for unknown words.
#'
#' Words that are flagged as \code{known} in the input will be associated with
#' \code{NA} in the output.
#'
#' Children acquire words that are hubs in their learning environment,
#' regardless of what they currently know [2].
#'
#' Growth value is the degree of the word in the presumed learning environment.
#'
#' @seealso \code{\link[netgrowr]{GrowthValue}}, \code{\link[netgrowr]{PreferentialAttachment}}, \code{\link[netgrowr]{LureOfTheAssociates}}, \code{\link[netgrowr]{opsahl_weighted_degree}}
#'
#' @references
#' 1. Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32(3), 245-251. doi:https://doi.org/10.1016/j.socnet.2010.03.006
#'
#' 2. Hills, T. T., Maouene, M., Maouene, J., Sheya, A., & Smith, L. (2009). Longitudinal analysis of early semantic networks: Preferential attachment or preferential acquisition? Psychological science, 20(6), 729-739. doi:10.1111/j.1467-9280.2009.02365.x
#'
#' @export
preferential_acquisition <- function(net, known, weighted = FALSE, alpha = 0.5) {
    pacq <- function(net, known) {
        gv <- colSums(net)
        gv[known] <- NA
        return(gv)
    }
    gv <- pacq(net > 0, known)
    if (weighted) {
        gvw <- pacq(net, known)
        gv <- opsahl_weighted_degree(gv, gvw, alpha)
    }
    names(gv) <- colnames(net)
    return(gv)
}

#' Reorder age of acquisition (AoA) values to align with network nodes
#'
#' @param net A binary adjacency matrix.
#' @param aoa Named vector of AoA values or a data frame with labels in column 1
#'   and AoA values in column 2.
#' @return A named vector reordered to align with network nodes (column names in
#'   adjacency matrix).
#'
align_aoa_vec <- function(net, aoa) {
    if (is.data.frame(aoa)) {
        aoa_vec <- aoa[[2]]
        names(aoa_vec) <- aoa[[1]]
    } else {
        aoa_vec <- aoa
    }
    return(aoa_vec[colnames(net)])
}

#' Re-weight degree as combination of weighted and unweighted values
#'
#' @param u Degree as derived from network with weights discarded (unweighted)
#' @param w Degree as derived from network with weights included (weighted)
#' @param alpha Mixing parameter, where zero means to only consider \code{u} and one means to only consider \code{w}. Use 0.5 to consult both sources equally.
#' @return Weighted degree
#'
#' @references
#' 1. Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32(3), 245-251. doi:https://doi.org/10.1016/j.socnet.2010.03.006
opsahl_weighted_degree <- function(u, w, alpha) {
    return(u^(1 - alpha) * w^alpha)
}

#' Flag words acquired before a given month
#'
#' @param month Positive integer indicating a reference month.
#' @param aoa Vector of age of acquisition values
#' @return Logical vector
#'
#' A word flagged as "known" entered the vocabulary prior to the indicated
#' month.
is_known <- function(month, aoa) {
    return(aoa < month)
}

#' Flag words acquired during a given month
#'
#' @param month Positive integer indicating a reference month.
#' @param aoa Vector of age of acquisition values
#' @return Logical vector
#'
#' A word flagged as "learned" enters the vocabulary in the indicated month.
is_learned <- function(month, aoa) {
    return(aoa == month)
}
