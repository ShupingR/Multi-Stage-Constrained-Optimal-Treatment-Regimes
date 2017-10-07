myp <- function (objectCDE, xNew, B = 1000, predictionBandProb = FALSE) 
{
  if (is.vector(xNew)) 
    xNew = as.matrix(xNew)
  if (class(objectCDE) != "FlexCoDE") 
    stop("Object should be of type FlexCoDE")
  zGrid = seq(from = 0, to = 1, length.out = B)
  if (is.null(objectCDE$bestI)) 
    objectCDE$bestI = objectCDE$nIMax
  coeff = predict(objectCDE$regressionObject, xNew, maxTerms = objectCDE$bestI)
  basisZNew = calculateBasis(zGrid, objectCDE$bestI, objectCDE$system)
  estimates = coeff %*% t(basisZNew)
  binSize = (1)/(B + 1)
  delta = ifelse(!is.null(objectCDE$bestDelta), objectCDE$bestDelta, 
                 0)
  estimates = t(apply(estimates, 1, function(xx) .normalizeDensity(binSize, 
                                                                   xx, delta)))
  estimates = estimates/(objectCDE$zMax - objectCDE$zMin)
  returnValue = NULL
  returnValue$CDE = estimates
  returnValue$z = seq(from = objectCDE$zMin, to = objectCDE$zMax, 
                      length.out = B)
  if (predictionBandProb == FALSE) 
    return(returnValue)
  th = matrix(NA, nrow(returnValue$CDE), 1)
  for (i in 1:nrow(returnValue$CDE)) {
    th[i] = .findThresholdHPD((objectCDE$zMax - objectCDE$zMin)/B, 
                              returnValue$CDE[i, ], predictionBandProb)
  }
  returnValue$th = th
  return(returnValue)
}