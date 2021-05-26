#' Reassign cells after running slingshot with reassign = FALSE
#'
#' @description Will return the weights after a final reassignment
#' @param sds Slingshot datasets
#' @return the weights
#' @importFrom slingshot slingCurveWeights
#' @importFrom matrixStats rowMaxs rowMins
#' @export
sling_reassign <- function(sds) {
  # from slingshot package
  W <- slingCurveWeights(sds)
  D <- vapply(slingCurves(sds), function(p){ p$dist_ind }, rep(0, nrow(sds)))
  ordD <- order(D)
  W.prob <- W / rowSums(W)
  WrnkD <- cumsum(W.prob[ordD]) / sum(W.prob)
  Z <- D
  Z[ordD] <- WrnkD
  Z.prime <- 1 - Z^2
  Z.prime[W == 0] <- NA
  W0 <- W
  W <- Z.prime / matrixStats::rowMaxs(Z.prime,na.rm = TRUE) #rowMins(D) / D
  W[is.nan(W)] <- 1 # handle 0/0
  W[is.na(W)] <- 0
  W[W > 1] <- 1
  W[W < 0] <- 0
  W[W0 == 0] <- 0
  # add if z < .5
  idx <- Z < .5
  W[idx] <- 1 #(rowMins(D) / D)[idx]
  # drop if z > .9 and w < .1
  ridx <- rowMaxs(Z, na.rm = TRUE) > .9 &
    rowMins(W, na.rm = TRUE) < .1
  W0 <- W[ridx, ]
  Z0 <- Z[ridx, ]
  W0[!is.na(Z0) & Z0 > .9 & W0 < .1] <- 0
  W[ridx, ] <- W0
  return(W)
}
