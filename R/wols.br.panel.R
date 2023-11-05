###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Panel Data


wols.br.panel <- function(deltaY, D, int.cov, pscore, i.weights, nthreads){
  #-----------------------------------------------------------------------------

  i.weights <- as.vector(i.weights * pscore/(1 - pscore))

  #Run weighted OLS
  beta.cal <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
                            subset = D==0,
                            weights = i.weights))

  tmp <- data.frame(deltaY,int.cov, D)
  names(tmp)[which(names(tmp) == "X.Intercept.")] <- "Intercept" # name gets transfomred, but must be "intercept"
  form <- paste0("deltaY~", paste(colnames(int.cov), collapse = "+"), "-1") # render formula, intercept is included in int.cov
  beta.cal <- suppressWarnings(fixest::feols(as.formula(form), subset = D==0,
                                            weights = i.weights, data = tmp, lean = TRUE))

  beta.cal <- suppressWarnings(beta.cal$coefficients)
  beta.cal <- as.vector(beta.cal)

  if(anyNA(beta.cal)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }

  #get fitted values
  #out.delta <-  as.numeric(tcrossprod(beta.cal, int.cov))
  out.delta <- cpppar_xbeta(int.cov, matrix(beta.cal, nrow = 1, ncol = length(beta.cal)), nthreads = nthreads) # equivalent to tcrossprod(beta.cal, int.cov))
  # return fitted values
  return(list(out.reg = out.delta))

}
