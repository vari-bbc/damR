#' QRILC imputation
#'
#' @param filteredData Missingness filtered data
#' @param tune.sigma coefficient that controls the sd of the MNAR distribution
#    - tune.sigma = 1 if the complete data distribution is supposed to be gaussian
#    - 0 < tune.sigma < 1 if the complete data distribution is supposed to be left-censored
#'
#' @returns Matrix with NA values replaced by imputed values
#'
#' @importFrom collapse fgroup_by fmutate fnunique fungroup
#' @importFrom stats lm qnorm quantile
#'
#' @export
impute.QRILC <- function(filteredData, tune.sigma = 1) {
  value <- imputed <- Sample <- NULL
  fgroup_by(filteredData, "Metabolite") |>
    fmutate(
      "imputed" = .imputer(value, nFeatures = fnunique(Sample), tune.sigma = tune.sigma),
      "diff" = value - imputed,
      "value" = ifelse(is.na(value), imputed, value)
    ) |>
    fungroup()
}

.imputer <- function(metabolite, nFeatures, tune.sigma) {
  pNAs <- mean(is.na(metabolite))
  if (pNAs == 0) {
    return(metabolite)
  }
  upper.q <- 0.99
  q.normal <- qnorm(
    seq((pNAs + 0.001), (upper.q + 0.001), (upper.q - pNAs) / (upper.q * 100)),
    mean = 0,
    sd = 1
  )
  q.curr.sample <- quantile(
    metabolite,
    probs = seq(0.001, (upper.q + 0.001), 0.01),
    na.rm = TRUE
  )
  lm.qr <- lm(q.curr.sample ~ q.normal)
  mean.CDD <- lm.qr$coefficients[1]
  sd.CDD <- as.numeric(lm.qr$coefficients[2])
  tmvtnorm::rtmvnorm(
    n = nFeatures,
    mean = mean.CDD,
    sigma = sd.CDD * tune.sigma,
    upper = qnorm((pNAs + 0.001), mean = mean.CDD, sd = sd.CDD),
    algorithm = c("gibbs")
  )
}
