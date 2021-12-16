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
#' @rdname icp.torus
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

#' @param object \code{icp.torus} object
#' @param newdata \code{data}
#' @param ... additional parameter
#'
#' @rdname icp.torus
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


#' Plot icp.torus object
#'
#' \code{plot.icp.torus} plots \code{icp.torus} object with some options.
#'
#' @param x \code{icp.torus} object
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#'   or \eqn{[-\pi, \pi)^d}. Default is \code{NULL}.
#' @param level either a scalar or a vector, or even \code{NULL}. Default value
#'   is 0.1.
#' @param ellipse A boolean index which determines whether plotting ellipse-intersections. Default is \code{TRUE}.
#' @param out An option for returning the ggplot object. Default is \code{FALSE}.
#' @param type A string. One of "mix", "max" or "e". This argument is only available if \code{icp.torus}
#'   object is fitted with method = "mixture". Default is \code{NULL}. If \code{type != NULL}, argument
#'   \code{ellipse} automatically becomes \code{FALSE}. If "mix", it plots based on von Mises mixture.
#'   If "max", it plots based on von Mises max-mixture. If "e", it plots based on ellipse-approximation.
#' @param ... additional parameters. For plotting icp.torus, these parameters are for ggplot2::ggplot().
#'
#' @rdname icp.torus
#' @method plot icp.torus
#' @export
plot.icp.torus <- function(x, data = NULL, level = 0.1, ellipse = TRUE, out = FALSE, type = NULL, ...){

  icp.torus <- x
  if (!is.null(type) && !(type %in% c("mix", "max", "e"))) {stop("invalid input: type must be one of \"mix\", \"max\" or \"e\".")}
  if (!is.null(type)) {ellipse <- FALSE}
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
      if (type == "mix"){b <- data.frame(ia$eval.point, ia$Chat_mix == 1)}
      else if (type == "max"){b <- data.frame(ia$eval.point, ia$Chat_max == 1)}
      else {b <- data.frame(ia$eval.point, ia$Chat_e == 1)}
    } else {
      b <- data.frame(ia$eval.point, ia$Chat_kmeans == 1)
    }
    colnames(b) <- c("phi","psi", method)
    g0 <- ggplot2::ggplot(...) + ggplot2::geom_contour(ggplot2::aes(x = .data$phi, y = .data$psi, z = ifelse(b[, 3], 1, 0)),
                                                    data = b, size = 1, lineend = "round" ) +
      ggplot2::geom_point(mapping = ggplot2::aes(x = .data$x, y = .data$y), data = data.frame(x = data[, 1],y = data[, 2])) +
      ggplot2::scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi)) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi)) +
      ggplot2::ggtitle(paste0("(Inductive) conformal prediction set with level = ", (1-level) * 100,"%"))
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
# 1. print, 2. plot

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

#' plot the conformal prediction set based on kernel density estimation
#'
#' @param x \code{cp.torus.kde} object
#' @param level.id an integer among \code{1:length(cp.torus$level)}.
#' @param ... additional parameter for ggplot2::ggplot()
#'
#' @rdname cp.torus.kde
#' @method plot cp.torus.kde
#' @export

plot.cp.torus.kde <- function(x, level.id = 1, ...){
  #  level.id is an integer among 1:length(cp.torus$level).

  if (is.null(x$cp.torus)) {stop("Conformal prediction set has not been evaluated.")}

  cp.torus <- x$cp.torus[x$cp.torus$level == x$level[level.id], ]

  g.new = ggplot2::ggplot(...) +
    ggplot2::geom_contour(ggplot2::aes(.data$phi, .data$psi, z = ifelse(.data$Cn, 1, 0)), data = cp.torus,
                          size = 1, lineend = "round" ) +
    ggplot2::geom_point(mapping = ggplot2::aes(.data$phi, .data$psi), data = as.data.frame(x$data.sorted)) +
    ggplot2::scale_x_continuous(breaks = c(0, 1, 2, 3, 4) * pi / 2,
                                labels = c("0","pi/2", "pi", "3pi/2", "2pi"), limits = c(0, 2 * pi)) +
    ggplot2::scale_y_continuous(breaks = c(0, 1, 2, 3, 4) * pi / 2,
                                labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi)) +
    ggplot2::ggtitle(label = paste0("Conformal prediction set with level = ", (1 - x$level[level.id]) * 100, "%"))
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


#' Plot cluster.obj object
#'
#' \code{plot.clus.torus} plots clustering results, which is given by \code{cluster.obj} object, with some options.
#'
#' @param x \code{cluster.obj} object
#' @param assignment A string. One of "outlier", "log.density", "posterior", "mahalanobis". Default is "outlier".
#' @param out An option for returning the ggplot object. Default is \code{FALSE}.
#' @param overlay A boolean index which determines whether plotting ellipse-intersections on clustering plots.
#'   Default is \code{FALSE}.
#' @param ... additional parameter for ggplot2::ggplot()
#'
#' @rdname cluster.assign.torus
#' @method plot cluster.obj
#' @export
plot.cluster.obj <- function(x, assignment = "outlier", overlay = FALSE, out = FALSE, ...){

  cluster.obj <- x
  if(is.null(data)) {stop("invalid input: data must be input.")}
  if(!is.character(assignment)) {stop("invalid input: assignment must be a character, one of \"log.density\", \"posterior\", \"outlier\", \"mahalanobis\".")}
  if(is.null(assignment)) {assignment <- "outlier"}

  data <- on.torus(cluster.obj$data)
  d <- ncol(data)
  if (d != cluster.obj$icp.torus$d) {stop("dimension mismatch: dimension of icp.torus object and input data are not the same.")}
  if (assignment == "outlier") {membership <- cluster.obj$cluster.id.outlier}
  else if (assignment == "log.density") {membership <- cluster.obj$cluster.id.log.density}
  else if (assignment == "posterior") {membership <- cluster.obj$cluster.id.by.posterior}
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
    plot_list[[i]] <- ggplot2::ggplot(...) +
      ggplot2::scale_color_manual(values = c(grDevices::hcl(h = hues, l = 65, c = 100), "#999999")[1:colnum]) +
      ggplot2::geom_point(ggplot2::aes(x = .data$angle1, y = .data$angle2, color = .data$membership), data = plotdata) +
      ggplot2::scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi")) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi")) +
      ggplot2::coord_cartesian(xlim = c(0, 2*pi), ylim = c(0, 2*pi), expand = FALSE)
  }
  legend <- cowplot::get_legend(plot_list[[1]] + ggplot2::theme(legend.box.margin = ggplot2::margin(0,0,0,12)))

  # plotlist <- list()
  # for (i in 1:(d-1)){
  #   for (j in 1:(d-1)){
  #     if (j < i) {plotlist[[(d-1) * (i-1) + j]] <- NULL}
  #     else {plotlist[[(d-1) * (i-1) + j]] <- plot_list[[(i-1) * (2*(d-1)-i)/2 + j]] + ggplot2::theme(legend.position = "none") +
  #       ggplot2::xlab(label = NULL) + ggplot2::ylab(label = NULL)}
  #     if (i == 1) {plotlist[[(d-1) * (i-1) + j]] <- plotlist[[(d-1) * (i-1) + j]] + ggplot2::ylab(label = angle.names[j+1])}
  #     if (j == (d-1)) {plotlist[[(d-1) * (i-1) + j]] <- plotlist[[(d-1) * (i-1) + j]] + ggplot2::xlab(label = angle.names[i])}
  #   }
  # }
  if (overlay == TRUE){
    level <- cluster.obj$level
    n2 <- cluster.obj$icp.torus$n2
    ialpha <- ifelse((n2 + 1) * level < 1, 1, floor((n2 + 1) * level))
    t <- cluster.obj$icp.torus$score_ellipse[ialpha]
    ellipse.param <- cluster.obj$icp.torus$ellipsefit

    J <- length(ellipse.param$c)
    theta <- seq(0, 2 * pi, length.out = 199)
    Z <- cbind(cos(theta), sin(theta))

    shift <- matrix(0, ncol = 2, nrow = 9)
    shift[, 1] <- c(0, 2 * pi, -2 * pi)
    shift[, 2] <- rep(c(0, 2 * pi, -2 * pi), each = 3)

    for (i in 1:nrow(coord)){
      i <- i
      ang1 <- coord[i, 1]
      ang2 <- coord[i, 2]

      for (j in 1:J){
        mu <- ellipse.param$mu[j, c(ang1, ang2)]
        Sinv <- ellipse.param$Sigmainv[[j]][c(ang1, ang2), c(ang1, ang2)]
        c.minus.t <- ellipse.param$c[j] - t

        if (c.minus.t <= 0) {next}

        M <- eigen(Sinv/c.minus.t)
        Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
        R <- Mmhalf %*% t(Z)
        for( shift.id in 1:9){
          RR <- R + mu + shift[shift.id,]
          plot.data <- data.frame(angle1 = RR[1,], angle2 = RR[2,], value = 1)
          plot_list[[i]] <- plot_list[[i]] + ggplot2::geom_polygon(ggplot2::aes(x = .data$angle1, y = .data$angle2),
                                                                   color = "blue",alpha = 0.1,
                                                                   data = plot.data)
        }
      }
      # plot_list[[i]] <- plot_list[[i]] + local({
      #   i <- i
      #   ang1 <- coord[i, 1]
      #   ang2 <- coord[i, 2]
      #
      #   for (j in 1:J){
      #     mu <- ellipse.param$mu[j, c(ang1, ang2)]
      #     Sinv <- ellipse.param$Sigmainv[[j]][c(ang1, ang2), c(ang1, ang2)]
      #     c.minus.t <- ellipse.param$c[j] - t
      #
      #     if (c.minus.t <= 0) {next}
      #
      #     M <- eigen(Sinv/c.minus.t)
      #     Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
      #     R <- Mmhalf %*% t(Z)
      #     for( shift.id in 1:9){
      #       RR <- R + mu + shift[shift.id,]
      #       plot.data <- data.frame(angle1 = RR[1,], angle2 = RR[2,], value = 1)
      #       g2 <- ggplot2::geom_polygon(ggplot2::aes(x = .data$angle1, y = .data$angle2),
      #                                        color = "blue",alpha = 0.1,
      #                                        data = plot.data)
      #     }
      #   }
      #   g2 <- g2 + ggplot2::coord_cartesian(xlim = c(0, 2*pi), ylim = c(0, 2*pi), expand = FALSE)
      # })
    }
  }

  if (out == FALSE) {
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
    print(cowplot::plot_grid((cowplot::plot_grid(plotlist = plotlist, ncol = (d-1), byrow = FALSE)),
                                              legend, ncol = 2, rel_widths = c(1, 0.3)))
    }
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

#' @param object \code{kmeans.torus} object
#' @param newdata \code{data} to be input
#' @param ... additional parameter
#'
#' @rdname kmeans.torus
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

#' @param x \code{hyperparam.J} object
#' @param ... additional parameter for ggplot2::ggplot()
#'
#' @method plot hyperparam.J
#' @rdname hyperparam.J
#' @export
plot.hyperparam.J <- function(x, ...){

  hyperparam.J <- x
  ggplot2::ggplot(ggplot2::aes(x = .data$J, y = .data$criterion),
                  data = as.data.frame(hyperparam.J$IC.results), ...) + ggplot2::geom_point(alpha = 0.5, shape = 4) +
    ggplot2::geom_point(ggplot2::aes(x = .data$J, y = .data$criterion),
                        data = as.data.frame(hyperparam.J$IC.results[which.min(hyperparam.J$IC.results[, 2]), ]),
                        shape = 19, size = 1.5) +
    ggplot2::geom_line(alpha = 0.5) + ggplot2::ggtitle(paste0("criterion: ", hyperparam.J$criterion, ', chosen J=', hyperparam.J$Jhat))
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

#' @param x \code{hyperparam.alpha} object
#' @param ... additional parameter for ggplot2::ggplot()
#' @method plot hyperparam.alpha
#' @rdname hyperparam.alpha
#' @export
plot.hyperparam.alpha <- function(x, ...){

  hyperparam.alpha <- x
  nclusters.length <- rle(hyperparam.alpha$alpha.results$ncluster)$lengths
  length <- max(nclusters.length)
  length.index <- which.max(nclusters.length)
  length.sum <- ifelse(length.index == 1, 0, sum(nclusters.length[1:(length.index - 1)]))

  ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$ncluster),
                  data = as.data.frame(hyperparam.alpha$alpha.results), ...) + ggplot2::geom_point() +
    ggplot2::geom_point(ggplot2::aes(x = .data$alpha, y = .data$ncluster), data = hyperparam.alpha$alpha.results[(length.sum + 1):(length.sum + length), ],
                        color = "blue") +
    ggplot2::ggtitle(paste0("number of clusters, chosen alpha=", hyperparam.alpha$alphathat))
}

# hyperparam.torus object ------------------------------------------
#' @method print hyperparam.torus
#' @export
print.hyperparam.torus <- function(x, ...){

  hyperparam.torus <- x
  cat("Type of conformity score:", hyperparam.torus$method, "\n")
  cat("Optimizing method:", hyperparam.torus$option, "\n")
  cat("-------------\n")
  if (hyperparam.torus$method[1] == "kde"){
    cat("Optimally chosen parameters. Concentration = ", hyperparam.torus$khat,
        ", alpha = ", hyperparam.torus$alphahat, "\n")
  }else{
    cat("Optimally chosen parameters. Number of components = ", hyperparam.torus$Jhat,
        ", alpha = ", hyperparam.torus$alphahat, "\n")
  }
  # cat("Optimally chosen parameters",
  #     paste0("(", ifelse(hyperparam.torus$method[1] == "kde", "concentration", "J"),
  #     ", alpha):"),
  #     hyperparam.torus$optim$hyperparam, "\n")
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
#' @param x \code{hyperparam.torus} object
#' @param color A string for plotting \code{hyperparam.torus} object, whose criterion option is \code{option = "elbow"}.
#'   One of "auto", "sequential", or "qualitative". If \code{color = "auto"},
#'   color assignment will be done automatically based on the number of J or concentration.
#'   If \code{color = "sequential"}, color assignment will be done by regarding each J or concentration as quantitative variable.
#'   If \code{color = "qualitative"}, color assignment will be done by regarding each J or concentration as qualitative variable.
#'   Default is \code{color = "auto"}.
#' @param ... additional parameter for ggplot2::ggplot()
#' @rdname hyperparam.torus
#' @export
plot.hyperparam.torus <- function(x, color = "auto", ...){

  hyperparam.torus <- x
  if (hyperparam.torus$option == "elbow" && hyperparam.torus$method[1] == "kde"){
    data <- as.data.frame(hyperparam.torus$results)
    if ((color == "qualitative") || (color == "auto" && length(unique(hyperparam.torus$results)) < 10)){
      data$k <- as.factor(data$k)
    }
    ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$criterion, color = .data$k, type = .data$k),
                    data = data, ...) + ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(x = .data$alpha, y = .data$criterion, group = .data$k))
  } else if (hyperparam.torus$option == "elbow" && hyperparam.torus$method[1] != "kde"){
    data <- as.data.frame(hyperparam.torus$results)
    if ((color == "qualitative") || (color == "auto" && length(unique(hyperparam.torus$results)) < 10)){
      data$k <- as.factor(data$J)
    }
    ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$criterion, color = .data$J, type = .data$J),
                    data = data, ...) + ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(x = .data$alpha, y = .data$criterion, group = .data$J))
  } else if (hyperparam.torus$option != "elbow"){
    g1 <- ggplot2::ggplot(ggplot2::aes(x = .data$J, y = .data$criterion),
                          data = as.data.frame(hyperparam.torus$IC.results), ...) + ggplot2::geom_point(alpha = 0.5, shape = 4) +
      ggplot2::geom_point(ggplot2::aes(x = .data$J, y = .data$criterion),
                          data = as.data.frame(hyperparam.torus$IC.results[which.min(hyperparam.torus$IC.results[, 2]), ]),
                          shape = 19, size = 1.5) +
      ggplot2::geom_line(alpha = 0.5) +
      ggplot2::ggtitle(paste0("criterion: ", hyperparam.torus$option, ', chosen J=', hyperparam.torus$Jhat))

    nclusters.length <- rle(hyperparam.torus$alpha.results$ncluster)$lengths
    length <- max(nclusters.length)
    length.index <- which.max(nclusters.length)
    length.sum <- ifelse(length.index == 1, 0, sum(nclusters.length[1:(length.index - 1)]))

    g2 <- ggplot2::ggplot(ggplot2::aes(x = .data$alpha, y = .data$ncluster),
                          data = as.data.frame(hyperparam.torus$alpha.results), ...) + ggplot2::geom_point() +
      ggplot2::geom_point(ggplot2::aes(x = .data$alpha, y = .data$ncluster), data = hyperparam.torus$alpha.results[(length.sum + 1):(length.sum + length), ],
                          color = "blue") +
      ggplot2::ggtitle(paste0("number of clusters, chosen alpha=", hyperparam.torus$alphathat))
    cowplot::plot_grid(g1, g2, nrow = 1)
  }
}

# clus.torus object ------------------------------------------
#' @method print clus.torus
#' @export
print.clus.torus <- function(x, ...){
  print(x$cluster.obj)
}

#' Plot clus.torus object
#'
#' \code{plot.clus.torus} plots \code{clus.torus} object with some options.
#'
#' @param x \code{clus.torus} object
#' @param panel One of 1 or 2 which determines the type of plot. If \code{panel = 1},
#'   \code{x$cluster.obj} will be plotted, if \code{panel = 2}, \code{x$icp.torus} will be plotted.
#'   If \code{panel = 3}, \code{x$hyperparam.select} will be plotted. Default is \code{panel = 1}.
#' @param assignment A string. One of "outlier", "log.density", "posterior", "mahalanobis". Default is "outlier".
#' @param overlay A boolean index which determines whether plotting ellipse-intersections on clustering plots. Default is \code{FALSE}.
#'   Only available for \code{panel = 1}.
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#'   or \eqn{[-\pi, \pi)^d}. Default is \code{NULL}.
#' @param type A string. One of "mix", "max" or "e". This argument is only available if \code{icp.torus}
#'   object is fitted with method = "mixture". Default is \code{NULL}. If \code{type != NULL}, argument
#'   \code{ellipse} automatically becomes \code{FALSE}. If "mix", it plots based on von Mises mixture.
#'   If "max", it plots based on von Mises max-mixture. If "e", it plots based on ellipse-approximation.
#' @param ellipse A boolean index which determines whether plotting ellipse-intersections. Default is \code{TRUE}. Only available
#'   for \code{panel = 2}.
#' @param out An option for returning the ggplot object. Default is \code{FALSE}.
#' @param ... additional parameter for ggplot2::ggplot()
#' @method plot clus.torus
#' @rdname clus.torus
#' @export
plot.clus.torus <- function(x, panel = 1, assignment = "outlier", data = NULL,
                            ellipse = TRUE, type = NULL, overlay = FALSE, out = FALSE, ...){
  if (!(panel %in% c(1, 2, 3))) {stop("invalid input: panel must be one of 1, 2 or 3.")}
  if (panel == 1){plot(x$cluster.obj, assignment = assignment, overlay = overlay, out = out, ...)}
  else if (panel == 2){plot(x$icp.torus, data = data, level = x$cluster.obj$level, ellipse = ellipse, out = out, type = type, ...)}
  else if (panel == 3){plot(x$hyperparam.select, ...)}
}

