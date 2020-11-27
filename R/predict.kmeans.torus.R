
predict.kmeans.torus <- function(data, kmeans){
  extrinsic.results <- kmeans$extrinsic.results
  extrinsic.data <- cbind(cos(data), sin(data))

  predict.kmeans <- apply(extrinsic.data, 1, function(r)
    {which.min(colSums((t(extrinsic.results$centers) - r)^2))})

  return(predict.kmeans)
}
