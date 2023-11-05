###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


pscore.cal <- function(D, int.cov, i.weights, n, nthreads){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  tmp <- data.frame(D,int.cov)
  names(tmp)[which(names(tmp) == "X.Intercept.")] <- "Intercept" # name gets transfomred, but must be "intercept"
  form <- paste0("D~", paste(colnames(int.cov), collapse = "+"), "-1") # render formula, intercept is included in int.cov
  pslogit <- suppressWarnings(fixest::feglm(as.formula(form), family = "binomial",
                                         weights = i.weights, data = tmp, lean = TRUE))

  init.gamma <- suppressWarnings(pslogit$coefficients)
  init.gamma <- as.vector(init.gamma)

  rm(tmp)
  gc()

  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))

  flag <- ifelse(pscore.cal$converged, 0, 1)

  gamma.cal <- try(pscore.cal$argument)

  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {

    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)

    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }

  #Compute fitted pscore and weights for regression
  pscore.index <- cpppar_xbeta(int.cov, matrix(gamma.cal, nrow = 1, ncol = length(gamma.cal)), nthreads = nthreads) # equivalent to tcrossprod(gamma.cal, int.cov)
  pscore.index <- as.vector(pscore.index)
  pscore <- as.numeric(stats::plogis(pscore.index))

  if(flag==1) {
    warning("trust algorithm did not converge when estimating propensity score. \n Used IPT algorithm a la Graham et al (2012)")
  }
  if(flag==2) {
    warning(" Used glm algorithm to estimate propensity score as trust and IPT method did not converge")

    if(pslogit$converged == FALSE){
      warning(" glm algorithm did not converge")
    }

  }

  if(anyNA(pscore)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }

  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))

}
