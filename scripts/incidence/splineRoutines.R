
# t, time points for which to infer the derivative
# gamToDerive, the spline/gam for which the derivatives are inferred
#gamDerivative <- function(tt.df, gamToDerive, outputpath) {
computeSplineDerivativeTable <- function(t, gamToDerive) {
  library(mgcv)
  tt.df = data.frame(t=t)

  ## finite difference interval, delta t
  eps <- 1e-7

  #When type="lpmatrix" then a matrix is returned which yields the values of the linear predictor (minus any offset)
  #when postmultiplied by the parameter vector (in this case se.fit is ignored). The latter option is most useful for
  #getting variance estimates for quantities derived from the model: for example integrated quantities, or derivatives of smooths.
  X0 <- predict(gamToDerive, tt.df)
  #print(X0)
  X0_e <- exp(X0)
  #print(X0_e)
  # super mall interval for derivation
  tt_eps.df <- tt.df + eps

  X1 <- predict(gamToDerive, tt_eps.df)
  X1_e <- exp(X1)
  # finite difference approximation of first derivative
  # the design matrix
  Xp <- (X1_e - X0_e) / eps

  # first derivative
  d1_gam <- Xp#%*% coef(gamToDerive)

  d1_gam.table <- data.frame(t=t, value=exp(d1_gam))

  return(d1_gam.table)
}


#generate random values from a multivariate normal distribution
random_mvn <- function(n, mu, sig) {
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

# infer the confidence interval by taking the simultaneous interval for a function f(x) at a set of M locations
# simulation-based approach by Rupert et. al (2003)
calculateCI <- function(gam_mod_cs, xx) {
  #  Bayesian covariance matrix of the model coefficients Vb, with unknown (to be estimated) smoothing parameter
  Vb <- vcov(gam_mod_cs, unconditional = TRUE)

  # generate fitted values and CI in link function...
  splinePred_gam_cs <- predict(gam_mod_cs, data.frame(t=xx), se.fit = TRUE, type = "link")

  N <- 10000

  # N draws from approximately multivariate normaldistribution with mean vector 0 and covariance matrix Vb
  BUdiff <- random_mvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

  # basis function at x
  Cx <- predict(gam_mod_cs, data.frame(t=xx), type = "lpmatrix")
  #deviations between the fitted and true parameters
  simDev <- Cx %*% t(BUdiff)
  #find the absolute values of the standardized deviations from the true model
  absDev <- abs(sweep(simDev, 1, splinePred_gam_cs$se.fit, FUN = "/"))
  #maximum of the absolute standardized deviations at the grid of x values for each simulation
  max_abs_sd <- apply(absDev, 2L, max)
  #calculate the critical value for a 95% simultaneous confidence interval/band
  crit <- quantile(max_abs_sd, prob = 0.95, type = 8)

  # ...and backtransform onto response scale
  return(data.frame(t=xx, value=gam_mod_cs$family$linkinv(splinePred_gam_cs$fit), se=splinePred_gam_cs$se.fit,
                    lower = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit-2*splinePred_gam_cs$se.fit)),
                    upper = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit+2*splinePred_gam_cs$se.fit)),
                    lowerSim = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit-crit*splinePred_gam_cs$se.fit)),
                    upperSim = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit+crit*splinePred_gam_cs$se.fit))
  ))
}

computeSpline <- function(input.table) {
  library(mgcv)
  #input table must include a column with t, value and variance

  #number of data points
  N<-nrow(input.table)
  #weight: choose 1/variance, maybe better standard deviation?
  weights <- 1/(input.table$variance)
  #normalize splines by sum of weights
  weights = weights/sum(weights)

  # optimze k! (default k=10)
  # increase k=degree of freedom (number of base functions)
  # until edf doesn't change substantially
  k=5
  edf_act=1
  repeat {
    edf = edf_act
    k<-k+5
    #compute splines, s=smooth function
    gam_mod_cs <- gam(value ~ s(t,bs="cs",k=k),
                      weights = weights,
                      data=input.table, method="REML", family=gaussian(link="log"))
    edf_act<- summary(gam_mod_cs)$edf
    if((edf_act-edf)<0.5 || k>= nrow(input.table))
      break
  }

  #return spline model
  return(gam_mod_cs)
}

computeSplineTable <- function(input.table) {
  gam_mod_cs <- computeSpline(input.table)
  # vector for which the outcomes should be predicted
  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)

  gam.table<-calculateCI(gam_mod_cs, xx)
  return(gam.table)
}

addSplineValuesForTrueN <- function(input.table, gam.table) {
  # optimize k (degree of freedom)
  k=5
  edf_act=1
  repeat {
    edf = edf_act
    k<-k+5
    #trueN comes from poisson distribution, hence no negative values (family=poisson -> link=logs)
    trueN_gam_cs <- gam(round(trueN) ~ s(t,bs="cs", k=k),
                        data=input.table, method="REML", family = poisson())
    edf_act<- summary(trueN_gam_cs)$edf
    if((edf_act-edf)<0.5 || k>= min(30,nrow(input.table)))
      break
  }

  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)

  gam.table_true<-calculateCI(trueN_gam_cs, xx)
  gam.table["value_trueN"] <-gam.table_true$value
  gam.table["value_trueN_lower"] <- gam.table_true$lowerSe
  gam.table["value_trueN_upper"] <- gam.table_true$upperSe
  gam.table["value_trueN_lowerSim"] <- gam.table_true$lowerSim
  gam.table["value_trueN_upperSim"] <- gam.table_true$upperSim

  return(gam.table)
}

computeSplineNewCasesTable <- function(input.table) {
  library(mgcv)
  # optimize k (degree of freedom)
  k=5
  edf_act=1
  repeat {
    edf = edf_act
    k<-k+5
    #trueN comes from poisson distribution, hence no negative values (family=poisson -> link=logs)
    trueN_gam_cs <- gam(round(new_cases) ~ s(t,bs="cs", k=k),
                        data=input.table, method="REML", family = poisson())
    edf_act<- summary(trueN_gam_cs)$edf
    if((edf_act-edf)<0.5 || k>= min(30,nrow(input.table)))
      break
  }

  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)

  gam.table<-calculateCI(trueN_gam_cs, xx)

  return(gam.table)
}


computeInterpolation <- function (input.table, ts, weights=NULL) {
  N_sampling <- 1000
  p_sample <- 0.5
  eps <- 10e-10
  #smoothing window of one week
  width <- 7

  if(is.null(weights))
    weights=rep(1, nrow(input.table))

  samplings <- matrix(ncol=length(ts), nrow=0)
  devations <- matrix(ncol=length(ts), nrow=0)
  rzeros <- matrix(ncol=length(ts), nrow=0)
  #alternative
  samplings_minusEpsi <- matrix(ncol=length(ts), nrow=0)
  samplings_plusEpsi <- matrix(ncol=length(ts), nrow=0)
  for(i in seq(N_sampling)) {
    # sample for each row if it is in the sampling set
    #sub_input.table <- input.table[runif(nrow(input.table))<p_sample, ]
    sub_value.table <- input.table[sample(seq(nrow(input.table)), round(nrow(input.table)*p_sample), replace = F, prob=NULL),]
    sub_value.table <- sub_value.table[order(sub_value.table$t),]

    #### first interpolate, collect and smooth median afterwards
    interpol <- approx(sub_value.table$t, sub_value.table$value, xout=ts)
    samplings <-rbind(samplings, interpol$y)
    #samplings.df <- rbind(samplings.df, data.frame(t=interpol$x, value=interpol$y, simu=i))

    #compute derivative
    interpol_smoothed_minusEpsilon<-approx(sub_value.table$t, sub_value.table$value, xout=ts)
    interpol_smoothed_plusEpsilon<-approx(sub_value.table$t, sub_value.table$value, xout=ts+eps)
    # samplings_minusEpsi <- rbind(samplings_minusEpsi, interpol_smoothed_minusEpsilon$y)
    # samplings_plusEpsi <- rbind(samplings_plusEpsi, interpol_smoothed_plusEpsilon$y)
  }

  # taking the quantiles of all samplings
  quantiles <- apply(samplings, 2, quantile, c(0.05, 0.5, 0.95), na.rm=T)
  smoothed_quantiles <- t(apply(quantiles, 1, filter, filter=rep(1,width)/width, sides=2))
  df <- data.frame(t=ts, smoothMedian=smoothed_quantiles[2,], smooth5=smoothed_quantiles[1,], smooth95=smoothed_quantiles[3,])
                   #qtMedian=smoothed_qts[2,], qt5=smoothed_qts[1,], qt95=smoothed_qts[3,])
  return(df)
}

computeSmoothedInterpolation <- function (input.table, ts, weights=NULL) {
  if(is.null(weights))
    weights=rep(1, nrow(input.table))

  N_sampling <- 1000
  p_sample <- 0.5
  eps <- 10e-10
  # average number of samples per day times 7 for a weekly average
  # (rounding to the nearest odd number)
  width <- floor((nrow(input.table)/length(ts)*7*p_sample)/2)*2+1

  if(is.null(weights))
    weights=rep(1, nrow(input.table))

  samplings_smoothed <- matrix(ncol=length(ts), nrow=0)
  devations <- matrix(ncol=length(ts), nrow=0)
  rzeros <- matrix(ncol=length(ts), nrow=0)
  for(i in seq(N_sampling)) {
    # sample for each row if it is in the sampling set
    #sub_input.table <- input.table[runif(nrow(input.table))<p_sample, ]
    sub_value.table <- input.table[sample(seq(nrow(input.table)), round(nrow(input.table)*p_sample), replace = F, prob=NULL),]
    sub_value.table <- sub_value.table[order(sub_value.table$t),]
    #### smooth interpolation and take median afterwards
    # convolutional filter with windowsize width
    smoothedVal = as.vector(filter(sub_value.table$value, rep(1,width)/width, sides=2))
    smothedTime = as.vector(filter(sub_value.table$t, rep(1,width)/width, sides=2))

    #inearly interpolate given data points
    interpol_smoothed<-approx(smothedTime, smoothedVal, xout=ts)
    samplings_smoothed <- rbind(samplings_smoothed, interpol_smoothed$y)

    #compute derivative
    interpol_smoothed_minusEpsilon<-approx(smothedTime, smoothedVal, xout=ts)
    interpol_smoothed_PlusEpsilon<-approx(smothedTime, smoothedVal, xout=ts+eps)
    deriv <- (interpol_smoothed_PlusEpsilon$y - interpol_smoothed_minusEpsilon$y)/(eps)
    devations <- rbind(devations, deriv)
    # R0 estimate here - for each bootstrap
    rzero <- exp(deriv)
    rzeros <- rbind(rzeros, rzero)

    # qt=y(t)/y(t-1)

  }

  # taking the quantiles of all smoothed samplings
  quantiles_smoothed <- apply(samplings_smoothed, 2, quantile, c(0.05, 0.5, 0.95), na.rm=T)
  quantiles_derivation <- apply(devations, 2, quantile, c(0.05, 0.5, 0.95), na.rm=T)
  quantiles_rzeros <- apply(rzeros, 2, quantile, c(0.05, 0.5, 0.95), na.rm=T)
  df <- data.frame(t=ts, smoothMedian=quantiles_smoothed[2,], smooth5=quantiles_smoothed[1,], smooth95=quantiles_smoothed[3,],
             derivativeMedian=quantiles_derivation[2,], derivative5=quantiles_derivation[1,], derivative95=quantiles_derivation[3,],
             rzeroMedian=quantiles_rzeros[2,], rzero5=quantiles_rzeros[1,], rzero95=quantiles_rzeros[3,])
  return(df)
}

computeRatio <- function(values, values_trueN) {
  return(values/values_trueN)
}
