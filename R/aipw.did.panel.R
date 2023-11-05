###################################################################################
# DR DiD estimator for the ATT with panel Data


aipw.did.panel <- function(deltaY, D, ps, out.reg, i.weights){
  #-----------------------------------------------------------------------------
  # Compute the AIPW estimator
  w.treat <- i.weights * D
  w.cont <- i.weights * (1 - D) * ps / (1 - ps)

  aipw.1 <- collapse::fmean(w.treat * (deltaY - out.reg)) / collapse::fmean(w.treat)
  aipw.0 <- collapse::fmean(w.cont * (deltaY - out.reg)) / collapse::fmean(w.cont)

  aipw.att <- aipw.1 - aipw.0

  return(aipw.att)
}
