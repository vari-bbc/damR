#' QRILC imputation
#'
#' @param dataSet.mvs Missingness filtered numeric matrix from `RemoveMissingness()`
#' @param tune.sigma coefficient that controls the sd of the MNAR distribution
#                        - tune.sigma = 1 if the complete data distribution is supposed
#                          to be gaussian
#                        - 0 < tune.sigma < 1 if the complete data distribution is
#                          supposed to be left-censored
#'
#' @returns Matrix with NA values replaced by imputed values
#'
#' @export
impute.QRILC <- function(filtered.mat, tune.sigma = 1) {
  nFeatures <- dim(filtered.mat)[1]
  nSamples <- dim(filtered.mat)[2]
  imputed.mat <- apply(filtered.mat, 2, .imputer, nFeatures, nSamples, tune.sigma)
  na.locations <- which(is.na(filtered.mat), arr.ind = TRUE)
  filtered.mat[na.locations] <- imputed.mat[na.locations]
  filtered.mat
}

.imputer <- function(metabolite, nFeatures, nSamples, tune.sigma) {
  pNAs <- mean(is.na(metabolite))
  if (pNAs == 0) {
    return(metabolite)
  }
  upper.q = 0.99
  q.normal = qnorm(
    seq((pNAs + 0.001), (upper.q + 0.001), (upper.q - pNAs) / (upper.q * 100)),
    mean = 0,
    sd = 1
  )
  q.curr.sample = quantile(
    metabolite,
    probs = seq(0.001, (upper.q + 0.001), 0.01),
    na.rm = TRUE
  )
  lm.qr = lm(q.curr.sample ~ q.normal)
  mean.CDD = lm.qr$coefficients[1]
  sd.CDD = as.numeric(lm.qr$coefficients[2])
  tmvtnorm::rtmvnorm(
    n = nFeatures,
    mean = mean.CDD,
    sigma = sd.CDD * tune.sigma,
    upper = qnorm((pNAs + 0.001), mean = mean.CDD, sd = sd.CDD),
    algorithm = c("gibbs")
  )
}
