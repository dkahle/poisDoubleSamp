#' Compute the score confidence interval
#'
#' Compute Rao's score confidence interval of a two-sample Poisson rate with misclassified data given fallible and infallible datasets.
#' 
#' @param data the vector of counts of the fallible data (z11, z12, z21, z22) followed by the infallible data (m011, m012, m021, m022, y01, y02)
#' @param N1 the opportunity size of group 1 for the fallible data
#' @param N2 the opportunity size of group 2 for the fallible data
#' @param N01 the opportunity size of group 1 for the infallible data
#' @param N02 the opportunity size of group 2 for the infallible data
#' @param conf.level confidence level of the interval
#' @return a named vector containing the lower and upper bounds of the confidence interval
#' @export scoreCI
#' @examples
#' \dontrun{
#' 
#' # small example
#' z11 <- 34; z12 <- 35; N1 <- 10; 
#' z21 <- 22; z22 <- 31; N2 <- 10;
#' m011 <- 9; m012 <- 1; y01 <- 3; N01 <- 3;
#' m021 <- 8; m022 <- 8; y02 <- 2; N02 <- 3;
#' data <- c(z11, z12, z21, z22, m011, m012, m021, m022, y01, y02)
#' 
#' waldCI(data, N1, N2, N01, N02) 
#' margMLECI(data, N1, N2, N01, N02)
#' profMLECI(data, N1, N2, N01, N02)
#' approxMargMLECI(data, N1, N2, N01, N02)
#' 
#' 
#' # big example :
#' z11 <- 477; z12 <- 1025; N1 <- 16186;
#' z21 <- 255; z22 <- 1450; N2 <- 18811;
#' m011 <- 38;  m012 <- 90; y01 <- 15; N01 <- 1500; 
#' m021 <- 41; m022 <- 200; y02 <-  9; N02 <- 2500;
#' data <- c(z11, z12, z21, z22, m011, m012, m021, m022, y01, y02)
#' 
#' waldCI(data, N1, N2, N01, N02) 
#' margMLECI(data, N1, N2, N01, N02)
#' profMLECI(data, N1, N2, N01, N02)
#' approxMargMLECI(data, N1, N2, N01, N02)
#'
#'
#' }
#'
scoreCI <- function(data, N1, N2, N01, N02, conf.level = .95){
  
  mle <- unname(fullMLE(data, N1, N2, N01, N02))
  se <- sePhi(mle, data, N1, N2, N01, N02)
  alpha <- 1-conf.level
  
  c(
    l = mle[1] - qnorm(1-alpha/2) * se,
    u = mle[1] + qnorm(1-alpha/2) * se
  )
}






