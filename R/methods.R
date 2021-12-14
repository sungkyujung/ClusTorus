# icp.torus object -------------------------------------------------
# 1. print, 2. summary, logLik


#' @method print icp.torus
#' @export
print.icp.torus <- function(x, ...){

  icp.torus <- x
  n <- 30
  cat("Data splited proportion", paste0(sum(icp.torus$split.id == 1), ":", sum(icp.torus$split.id != 1)), "\n")
  cat("\n")

  cat("Used method:", icp.torus$method, "\n")
  if(icp.torus$method == "kde"){
    cat("Concentration:", icp.torus$concentration, "\n")
    cat("-------------\n")
    cat("Conformity scores based on kde: \n")
    n <- min(icp.torus$n2, n)
    print(stats::quantile(icp.torus$score))
    cat("\n")
    print(utils::head(icp.torus$score, n = n))
    cat("\n")
    cat(icp.torus$n2-n, "conformity scores are omitted.")
  } else if (icp.torus$method == "mixture"){
    J <- length(icp.torus$ellipsefit$c)
    cat("Fitting method:", icp.torus$fittingmethod, "\n")
    cat("Number of mixture components:", J, "\n")
    cat("-------------\n")
    n <- min(icp.torus$n2, n)
    cat("Conformity scores based on von Mises mixture: \n")
    print(stats::quantile(icp.torus$score))
    cat("\n")
    print(utils::head(icp.torus$score, n = n))
    cat("\n")
    cat("Conformity scores based on max-mixture: \n")
    print(stats::quantile(icp.torus$score_max))
    cat("\n")
    print(utils::head(icp.torus$score_max, n = n))
    cat("\n")
    cat("Conformity scores based on ellipsoids: \n")
    print(stats::quantile(icp.torus$score_ellipse))
    cat("\n")
    print(utils::head(icp.torus$score_ellipse, n = n))
    cat("\n")
    cat(icp.torus$n2-n, "conformity scores are omitted.")
  } else {
    J <- length(icp.torus$ellipsefit$c)
    cat("Fitting method:", icp.torus$fittingmethod, "\n")
    cat("Number of mixture components:", J, "\n")
    cat("\n")
    if (!is.null(icp.torus$ellipsefit$singular)){
      cat("Number of singular components:", length(icp.torus$ellipsefit$singular), "\n")
      cat("Singular compnents J =", icp.torus$ellipsefit$singular, "\n")
    }
    cat("-------------\n")
    n <- min(icp.torus$n2, n)
    cat("Conformity scores based on ellpsoids: \n")
    print(stats::quantile(icp.torus$score_ellipse))
    cat("\n")
    print(utils::head(icp.torus$score_ellipse, n = n))
    cat("\n")
    cat(icp.torus$n2-n, "conformity scores are omitted.")
  }

  cat("\n\n")
  cat("Available components:\n")
  print(names(icp.torus))
}

#' @method logLik icp.torus
#' @export

logLik.icp.torus <- function(object, ...){

  icp.torus <- object
  method <- icp.torus$method
  j <- length(icp.torus$ellipsefit$c)
  d <- icp.torus$d
  if (method == "kde") {stop("method \'kde\' is not supported.")}
  else if (method == "mixture"){
    loglik <- icp.torus$fit$loglkhd.seq[length(icp.torus$fit$loglkhd.seq)]

    mixturefitmethod <- icp.torus$fittingmethod
    if (mixturefitmethod == "circular"){
      k <- d + 2
    } else if (mixturefitmethod == "axis-aligned"){
      k <- 2 * d + 1
    } else if (mixturefitmethod == "general"){
      k <- (d + 1)*(d + 2)/2
    }

    df <- k * j - 1
  } else {
    loglik <- icp.torus$ellipsefit$loglkhd

    kmeansfitmethod <- icp.torus$fittingmethod
    if (kmeansfitmethod == "homogeneous-circular"){
      k <- d + 1
    } else if (kmeansfitmethod == "heterogeneous-circular"){
      k <- d + 2
    } else {k <- (d + 1)*(d + 2)/2}

    df <- k * j - 1
  }
  cat("\'log Lik.\'", loglik, paste0("(df=", df,")"))
}

#' description
#'
#' more details
#'
#' @param object \code{icp.torus} object
#' @param newdata \code{data}
#' @param ... additional parameter
#'
#' @rdname predict.icp.torus
#' @method predict icp.torus
#' @export
predict.icp.torus <- function(object, newdata, ...){
  icp.torus <- object
  if (is.null(newdata)) {stop("newdata must be inpt for evaluating their conformity scores.")}
  newdata <- as.matrix(newdata)
  if (icp.torus$method == "kde"){
    phat <- kde.torus(icp.torus$X1, newdata, concentration = icp.torus$concentration)
    result <- sort(phat)
  } else if (icp.torus$method == "mixture"){
    result <- list(score = c(), score_max = c(), score_ellipse = c())
    result$score <-  sort(BAMBI::dvmsinmix(newdata, kappa1 = icp.torus$fit$parammat[2, ],
                                               kappa2 = icp.torus$fit$parammat[3, ],
                                               kappa3 = icp.torus$fit$parammat[4, ],
                                               mu1 = icp.torus$fit$parammat[5, ],
                                               mu2 = icp.torus$fit$parammat[6, ],
                                               pmix = icp.torus$fit$parammat[1, ], log = FALSE))
    result$score_max <- sort(do.call(pmax, as.data.frame(phat.eval(newdata, icp.torus$fit$parammat))))
    result$score_ellipse <- sort(do.call(pmax, as.data.frame(ehat.eval(newdata, icp.torus$ellipsefit))))
  } else {
    result <- sort(do.call(pmax, as.data.frame(ehat.eval(newdata, icp.torus$ellipsefit))))
  }

  return(result)
}


#' description
#'
#' more details
#'
#' @param x \code{x} object
#' @param data \code{data}
#' @param level a scalar
#' @param ellipse A boolean index.
#' @param out A boolean index.
#' @param ... additional parameter
#'
#' @rdname plot.icp.torus
#' @method plot icp.torus
#' @export
plot.icp.torus <- function(x, data = NULL, level = 0.1, ellipse = TRUE, out = FALSE, ...){

  icp.torus <- x
  
  #if(is.null(data)) {stop("invalid input: data must be input.")}
  
  if(is.null(data)) { data <- x$model}
  data <- on.torus(data)
  
  d <- ncol(data)
  if (d != icp.torus$d) {stop("dimension mismatch: dimension of icp.torus object and input data are not the same.")}
  n2 <- icp.torus$n2
  names <- colnames(data)
  angle.names <- colnames(data)


  method <- icp.torus$method
  if (method == "kde"){
   # if(ellipse == TRUE) {warning("method kde does not support the option ellipse. Automatically plot on the grid.")}
    ellipse <- FALSE
  }

  if (ellipse == FALSE){
    if (d > 2) {stop("dimension d > 2 cases does not support the option ellipse = FLASE")}
    ia <- icp.torus.eval(icp.torus, level = level, eval.point = grid.torus(d = d, grid.size = 100))
    if (method == "kde"){
      b <- data.frame(ia$eval.point, ia$Chat_kde == 1)
    } else if (method == "mixture"){
      b <- data.frame(ia$eval.point, ia$Chat_mix == 1)
    } else {
      b <- data.frame(ia$eval.point, ia$Chat_kmeans == 1)
    }
    colnames(b) <- c("phi","psi", method)
    g0 <- ggplot2::ggplot() + ggplot2::geom_contour(ggplot2::aes(x = .data$phi, y = .data$psi, z = ifelse(b[, 3], 1, 0)),
                                                    data = b, size = 1, lineend = "round" ) +
      ggplot2::geom_point(mapping = ggplot2::aes(x = .data$x, y = .data$y), data = data.frame(x = data[, 1],y = data[, 2])) +
      ggplot2::scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
      ggplot2::scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
    if(out == FALSE) {print(g0)}
    else {return(g0)}

  } else {
    ialpha <- ifelse((n2 + 1) * level < 1, 1, floor((n2 + 1) * level))
    t <- icp.torus$score_ellipse[ialpha]
    ellipse.param <- icp.torus$ellipsefit
    plot_list <- ploting.ellipsoids(data, ellipse.param, t, coord = t(utils::combn(d, 2)))
    plotlist <- list()
    for (i in 1:(d-1)){
      for (j in 1:(d-1)){
        if (j < i) {plotlist[[(d-1) * (i-1) + j]] <- NULL}
        else {plotlist[[(d-1) * (i-1) + j]] <- plot_list[[(i-1) * (2*(d-1)-i)/2 + j]] +
          ggplot2::xlab(label = NULL) + ggplot2::ylab(label = NULL)}
        if (i == 1) {plotlist[[(d-1) * (i-1) + j]] <- plotlist[[(d-1) * (i-1) + j]] + ggplot2::ylab(label = angle.names[j+1])}
        if (j == (d-1)) {plotlist[[(d-1) * (i-1) + j]] <- plotlist[[(d-1) * (i-1) + j]] + ggplot2::xlab(label = angle.names[i])}
      }
    }
    if (out == FALSE) {print(cowplot::plot_grid(plotlist = plotlist, ncol = (d-1), byrow = FALSE))}
    else {return(plot_list)}
  }
}
# icp.torus.eval object -------------------------------------------
# 1. print

#' @method print icp.torus.eval
#' @export
print.icp.torus.eval <- function(x, ...){

  icp.torus.eval <- x
  n = 10
  n <- min(nrow(icp.torus.eval$eval.point), n)
  if (!is.null(icp.torus.eval$Chat_kmeans)){
    cat("Conformal prediction set (Chat_kmeans)\n\n")
    cat("Testing inclusion to the conformal prediction set with level =", paste0(icp.torus.eval$level, ":\n"))
    cat("-------------\n")
    print(utils::head(data.frame(icp.torus.eval$eval.point, inclusion = icp.torus.eval$Chat_kmeans), n = n))
  } else if (!is.null(icp.torus.eval$Chat_kde)){
    cat("Conformal prediction set (Chat_kde)\n\n")
    cat("Testing inclusion to the conformal prediction set with level =", paste0(icp.torus.eval$level, ":\n"))
    cat("-------------\n")
    print(utils::head(data.frame(icp.torus.eval$eval.point, inclusion = icp.torus.eval$Chat_kde), n = n))
  } else {
    cat("Conformal prediction sets (Chat_mix, Chat_max, Chat_e)\n\n")
    cat("Testing inclusion to the conformal prediction set with level =", paste0(icp.torus.eval$level, ":\n"))
    cat("-------------\n")
    print(utils::head(data.frame(icp.torus.eval$eval.point, inclusion.mix = icp.torus.eval$Chat_mix,
                          inclusion.max = icp.torus.eval$Chat_max, inclusion.e = icp.torus.eval$Chat_e), n = n))
  }
  cat("\n", nrow(icp.torus.eval$eval.point) - n, "rows are omitted.")
}

# cp.torus.kde object ----------------------------------------------
# 1. print, 2. plot (to be added)

#' @method print cp.torus.kde
#' @export
print.cp.torus.kde <- function(x, ...){

  cp.torus.kde <- x
  n_score <- 30
  n_test <- 10
  cat("Conformal prediction sets (Lminus, Cn, Lplus) based on kde\n\n")

  cat(" with concentration: ", cp.torus.kde$concentration, "\n")
  cat("-------------\n")
  cat("Conformity scores based on kde: \n")
  print(stats::quantile(cp.torus.kde$phat.data))
  cat("\n")
  n_score <- min(length(cp.torus.kde$phat.data), n_score)
  print(utils::head(cp.torus.kde$phat.data, n = n_score))
  cat("\n")
  cat(length(cp.torus.kde$phat.data)-n_score, "conformity scores are omitted.")
  cat("\n\n")

  cat("Conformity scores of evaluation points, based on kde: \n")
  print(stats::quantile(cp.torus.kde$phat.grid))
  cat("\n")
  n_score <- min(length(cp.torus.kde$phat.grid), n_score)
  print(utils::head(cp.torus.kde$phat.grid, n = n_score))
  cat("\n")
  cat(length(cp.torus.kde$phat.grid)-n_score, "conformity scores are omitted.")
  cat("\n\n")

  cat("Testing inclusion to the conformal prediction set with level =", cp.torus.kde$level, ":\n")
  cat("-------------\n")
  n_test <- min(nrow(cp.torus.kde$cp.torus), n_test)
  print(utils::head(data.frame(cp.torus.kde$cp.torus), n = n_test))
  cat("\n", nrow(cp.torus.kde$cp.torus) - n_test, "rows are omitted.")
}

#' description
#'
#' more details
#'
#' @param x \code{cp.torus.kde} object
#' @param level.id a scalar
#' @param ... additional parameter
#'
#' @rdname plot.cp.torus.kde
#' @method plot cp.torus.kde
#' @export

plot.cp.torus.kde <- function(x,level.id = 1,...){ 
  #  level.id is an integer among 1:length(cp.torus$level).
  
  if (is.null(x$cp.torus)) {stop("Conformal prediction set has not been evaluated.")}
  
  cp.torus <- dplyr::filter(x$cp.torus, level == x$level[level.id])
  
  g.new = ggplot2::ggplot() +
    ggplot2::geom_contour(ggplot2::aes(phi,psi, z = ifelse(Cn, 1, 0)), data = cp.torus,
                          size = 1, lineend = "round" ) +
    ggplot2::geom_point(mapping = ggplot2::aes(phi,psi), data = as.data.frame(x$data.sorted)) +
    ggplot2::scale_x_continuous(breaks = c(0, 1, 2, 3, 4) * pi / 2,
                                labels = c("0","pi/2", "pi", "3pi/2", "2pi"), limits = c(0, 2 * pi)) +
    ggplot2::scale_y_continuous(breaks = c(0, 1, 2, 3, 4) * pi / 2,
                                labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
  g.new
}


# cluster.obj object -----------------------------------------------
# 1. print, 2. plot

#' @method print cluster.obj
#' @export
print.cluster.obj <- function(x, ...){

  cluster.obj <- x
  n <- 10
  n <- min(length(cluster.obj$cluster.id.log.density), n)
  cat("Number of clusters:", cluster.obj$ncluster, "\n")
  cat("-------------\n")
  cat("Clustering results by log density:\n")
  print(utils::head(cluster.obj$cluster.id.log.density, n = n))
  cat("cluster sizes:", purrr::map_int(1:cluster.obj$ncluster,
                                       function(x){sum(cluster.obj$cluster.id.log.density == x)}), "\n")
  cat("\n")
  cat("Clustering results by posterior:\n")
  print(utils::head(cluster.obj$cluster.id.by.posterior, n = n))
  cat("cluster sizes:", purrr::map_int(1:cluster.obj$ncluster,
                                       function(x){sum(cluster.obj$cluster.id.by.posterior == x)}), "\n")
  cat("\n")
  cat("Clustering results by representing outliers:\n")
  print(utils::head(cluster.obj$cluster.id.outlier, n = n))
  cat("cluster sizes:", purrr::map_int(1:(cluster.obj$ncluster+1),
                                       function(x){sum(cluster.obj$cluster.id.outlier == x)}), "\n")
  cat("\n Note: cluster id =", paste0(max(cluster.obj$cluster.id.outlier)), "represents outliers.\n")
  cat("\n")
  cat("Clustering results by Mahalanobis distance:\n")
  print(utils::head(cluster.obj$cluster.id.by.Mah.dist, n = n))
  cat("cluster sizes:", purrr::map_int(1:cluster.obj$ncluster,
                                       function(x){sum(cluster.obj$cluster.id.by.Mah.dist == x)}), "\n")
  cat("\n")
  cat(length(cluster.obj$cluster.id.log.density) - n, "clustering results are omitted.")
}


#' description
#'
#' more details
#'
#' @param x \code{x} object
#' @param method A string. One of "outlier", "log.density", "posterior", "mahalanobis". Default is "outlier".
#' @param out An option for returning the ggplot object.
#' @param ... additional parameter
#'
#' @rdname plot.cluster.obj
#' @method plot cluster.obj
#' @export
plot.cluster.obj <- function(x, method = c("outlier", "log.density", "posterior", "mahalanobis"), out = FALSE, ...){

  cluster.obj <- x
  if(is.null(data)) {stop("invalid input: data must be input.")}
  if(!is.character(method)) {stop("invalid input: method must be a character, one of \"log.density\", \"posterior\", \"outlier\", \"mahalanobis\".")}
  if(is.null(method)) {method <- "outlier"}

  method <- match.arg(method)
  data <- on.torus(cluster.obj$data)
  d <- ncol(data)
  if (d != cluster.obj$icp.torus$d) {stop("dimension mismatch: dimension of icp.torus object and input data are not the same.")}
  if (method == "outlier") {membership <- cluster.obj$cluster.id.outlier}
  else if (method == "log.density") {membership <- cluster.obj$cluster.id.log.density}
  else if (method == "posterior") {membership <- cluster.obj$cluster.id.by.posterior}
  else {membership <- cluster.obj$cluster.id.by.Mah.dist}

  coord <- t(utils::combn(d, 2))
  if (max(coord) > d | min(coord) < 1)
  {stop("invalid coordinates: out of dimensionality.")}

  outlier_id <- max(cluster.obj$cluster.id.outlier)

  unique.membership <- unique(membership)
  unique.membership <- sort(unique.membership[unique.membership != outlier_id])
  colnum <- length(membership)
  angle.names <- colnames(data)
  plot_list <- list()

  hues = seq(15, 375, length = outlier_id + 1)[unique.membership]
  for (i in 1:nrow(coord)){
    plotdata <- data.frame(angle1 = data[ , coord[i, 1]], angle2 = data[ , coord[i, 2]])
    membership[membership == outlier_id] <- "out"
    plotdata <- data.frame(plotdata, membership = factor(membership,
                                                         levels = c(as.character(unique.membership), "out")))
    plot_list[[i]] <- ggplot2::ggplot(ggplot2::aes(x = .data$angle1, y = .data$angle2, color = .data$membership), data = plotdata, ...) +
      ggplot2::scale_color_manual(values = c(grDevices::hcl(h = hues, l = 65, c = 100), "#999999")[1:colnum]) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi)) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
  }
  legend <- cowplot::get_legend(plot_list[[1]] + ggplot2::theme(legend.box.margin = ggplot2::margin(0,0,0,12)))

  plotlist <- list()
  for (i in 1:(d-1)){
    for (j in 1:(d-1)){
      if (j < i) {plotlist[[(d-1) * (i-1) + j]] <- NULL}
      else {plotlist[[(d-1) * (i-1) + j]] <- plot_list[[(i-1) * (2*(d-1)-i)/2 + j]] + ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(label = NULL) + ggplot2::ylab(label = NULL)}
      if (i == 1) {plotlist[[(d-1) * (i-1) + j]] <- plotlist[[(d-1) * (i-1) + j]] + ggplot2::ylab(label = angle.names[j+1])}
      if (j == (d-1)) {plotlist[[(d-1) * (i-1) + j]] <- plotlist[[(d-1) * (i-1) + j]] + ggplot2::xlab(label = angle.names[i])}
    }
  }
  if (out == FALSE) {print(cowplot::plot_grid((cowplot::plot_grid(plotlist = plotlist, ncol = (d-1), byrow = FALSE)),
                                              legend, ncol = 2, rel_widths = c(1, 0.3)))}
  else {return(plot_list)}
}

# kmeans.torus object ----------------------------------------------
# 1. print, 2. predict

#' @method print kmeans.torus
#' @export
print.kmeans.torus <- function(x, ...){

  kmeans.torus <- x
  n <- 30
  cat("Extrinsic k-means clustering with", length(kmeans.torus$size), "clusters of sizes", kmeans.torus$size, "\n")
  cat("\n")
  cat("Cluster means:\n")
  print(kmeans.torus$centers)
  cat("\n")
  cat("Clustering results:\n")
  n <- min(length(kmeans.torus$membership), n)
  print(utils::head(kmeans.torus$membership, n = n))
  cat(length(kmeans.torus$membership) - n, "clustering results are omitted.")
  cat("\n")
  cat("Within cluster sum of squares by cluster:\n")
  print(kmeans.torus$withinss)
  cat("\n")
  cat("Between cluster sum of squares by cluster:\n")
  print(kmeans.torus$betweenss)
}


#' description
#'
#' more details
#'
#' @param object \code{kmeans.torus} object
#' @param newdata \code{data} to be input
#' @param ... additional parameter
#'
#' @rdname predict.kmeans.torus
#' @method predict kmeans.torus
#' @export
predict.kmeans.torus <- function(object, newdata, ...){
  kmeans.torus <- object
  data <- on.torus(newdata)

  extrinsic.results <- kmeans.torus$extrinsic.results
  extrinsic.data <- cbind(cos(data), sin(data))

  pred.kmeans <- apply(extrinsic.data, 1, function(r)
  {which.min(colSums((t(extrinsic.results$centers) - r)^2))})

  return(pred.kmeans)
}

# hyperparam.J object ----------------------------------------------
# 1. print,  2. plot

#' @method print hyperparam.J
#' @export
print.hyperparam.J <- function(x, ...){
  hyperparam.J <- x
  cat("Results based on criterion", hyperparam.J$criterion, ":\n")
  print(hyperparam.J$IC.results)

  cat("Optimally chosen J:", hyperparam.J$Jhat, "\n")
  cat("\n")
  cat("Available components:\n")
  print(names(hyperparam.J))
}

#' @method plot hyperparam.J
#' @export
plot.hyperparam.J <- function(x, ...){

  hyperparam.J <- x
  ggplot2::ggplot(ggplot2::aes(x = .data$J, y = .data$criterion),
                  data = as.data.frame(hyperparam.J$IC.results)) + ggplot2::geom_point(alpha = 0.5, shape = 4) +
    ggplot2::geom_point(ggplot2::aes(x = .data$J, y = .data$criterion),
                        data = as.data.frame(hyperparam.J$IC.results[which.min(hyperparam.J$IC.results[, 2]), ]),
                        shape = 19, size = 1.5) +
    ggplot2::geom_line(alpha = 0.5) + ggplot2::ggtitle(paste("criterion: ", hyperparam.J$criterion))
}

# hyperparam.alpha object ------------------------------------------
# 1. print, 2. plot

#' @method print hyperparam.alpha
#' @export
print.hyperparam.alpha <- function(x, ...){
  hyperparam.alpha <- x
  cat("Clustering results by varying alpha on ", paste0("[", round(min(hyperparam.alpha$alpha.results$alpha), 5),
                                                        round(max(hyperparam.alpha$alpha.results$alpha), 5),"]"), ":\n")
  print(hyperparam.alpha$alpha.results)

  cat("Optimally chosen alpha:", hyperparam.alpha$alphahat, "\n")
}


#' @method plot hyperparam.alpha
#' @export
plot.hyperparam.alpha <- function(x, ...){

  hyperparam.alpha <- x
  nclusters.length <- rle(hyperparam.alpha$alpha.results$ncluster)$lengths
  length <- max(nclusters.length)
  length.index <- which.max(nclusters.length)
  length.sum <- ifelse(length.index == 1, 0, sum(nclusters.length[1:(length.index - 1)]))

  ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$ncluster),
                  data = as.data.frame(hyperparam.alpha$alpha.results)) + ggplot2::geom_point() +
    ggplot2::geom_point(ggplot2::aes(x = .data$alpha, y = .data$ncluster), data = hyperparam.alpha$alpha.results[(length.sum + 1):(length.sum + length), ],
                        color = "blue") +
    ggplot2::ggtitle("number of clusters")
}

# hyperparam.torus object ------------------------------------------
#' @method print hyperparam.torus
#' @export
print.hyperparam.torus <- function(x, ...){

  hyperparam.torus <- x
  cat("Type of conformity score:", hyperparam.torus$method, "\n")
  cat("Optimizing method:", hyperparam.torus$option, "\n")
  cat("-------------\n")
  cat("Optimally chosen parameters", paste0("(", ifelse(hyperparam.torus$method[1] == "kde", "concentration", "J"),
      ", alpha):"),
      hyperparam.torus$optim$hyperparam, "\n")
  if (hyperparam.torus$option == "elbow"){
    n <- min(15, nrow(hyperparam.torus$results))
    cat("Results based on criterion", hyperparam.torus$option, ":\n")
    print(utils::head(hyperparam.torus$results[sort(hyperparam.torus$results$criterion, index.return = TRUE)$ix, ], n))
    cat("\n")
    cat(nrow(hyperparam.torus$results) - n, "rows are omitted.\n")
  } else {
    cat("Results based on criterion", hyperparam.torus$option, ":\n")
    print(hyperparam.torus$IC.results)
    cat("\n")
    cat("Clustering results by varying alpha on", paste0("[", round(min(hyperparam.torus$alpha.results$alpha), 5),
                                                          round(max(hyperparam.torus$alpha.results$alpha), 5),"]"), ":\n")
    print(hyperparam.torus$alpha.results)
  }
  cat("\n")
  cat("Available components:\n")
  print(names(hyperparam.torus))
}

#' @method plot hyperparam.torus
#' @export
plot.hyperparam.torus <- function(x, ...){

  hyperparam.torus <- x
  if (hyperparam.torus$option == "elbow" && hyperparam.torus$method[1] == "kde"){
    ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$criterion, color = .data$k, type = .data$k),
                    data = as.data.frame(hyperparam.torus$results)) + ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(x = .data$alpha, y = .data$criterion, group = .data$k))
  } else if (hyperparam.torus$option == "elbow" && hyperparam.torus$method[1] != "kde"){
    ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$criterion, color = .data$J, type = .data$J),
                    data = as.data.frame(hyperparam.torus$results)) + ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(x = .data$alpha, y = .data$criterion, group = .data$J))
  } else if (hyperparam.torus$option != "elbow"){
    g1 <- ggplot2::ggplot(ggplot2::aes(x = .data$J, y = .data$criterion),
                          data = as.data.frame(hyperparam.torus$IC.results)) + ggplot2::geom_point(alpha = 0.5, shape = 4) +
      ggplot2::geom_point(ggplot2::aes(x = .data$J, y = .data$criterion),
                          data = as.data.frame(hyperparam.torus$IC.results[which.min(hyperparam.torus$IC.results[, 2]), ]),
                          shape = 19, size = 1.5) +
      ggplot2::geom_line(alpha = 0.5) + ggplot2::ggtitle(paste("option = ", hyperparam.torus$option))

    nclusters.length <- rle(hyperparam.torus$alpha.results$ncluster)$lengths
    length <- max(nclusters.length)
    length.index <- which.max(nclusters.length)
    length.sum <- ifelse(length.index == 1, 0, sum(nclusters.length[1:(length.index - 1)]))

    g2 <- ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$ncluster),
                          data = as.data.frame(hyperparam.torus$alpha.results)) + ggplot2::geom_point() +
      ggplot2::geom_point(ggplot2::aes(x = .data$alpha, y = .data$ncluster), data = hyperparam.torus$alpha.results[(length.sum + 1):(length.sum + length), ],
                          color = "blue") +
      ggplot2::ggtitle("number of clusters")
    cowplot::plot_grid(g1, g2, nrow = 1)
  }
}
