##### Those functions belong to the MSBVAR package
##### They are the property of Patrick Brandt

#a02mcmc
a02mcmc <- function (x) 
{
  return(mcmc(matrix(x$A0.posterior$A0, nrow = x$N2, ncol = length(x$A0.posterior$struct), 
                     byrow = T), thin = x$thin))
}

#cf.forecasts
cf.forecasts <- function (m1, m2) 
{
  rmse <- rmse(m1, m2)
  mae <- mae(m1, m2)
  return(c(rmse, mae))
}

# cf.forecasts
cf.forecasts <- function (m1, m2) 
{
  rmse <- rmse(m1, m2)
  mae <- mae(m1, m2)
  return(c(rmse, mae))
}

  
# ddirichlet
ddirichlet <- function (x, alpha) 
{
  dirichlet1 <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- sum((alpha - 1) * log(x))
    exp(sum(s) - logD)
  }
  if (!is.matrix(x)) 
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    else x <- t(x)
    if (!is.matrix(alpha)) 
      alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                      byrow = TRUE)
    if (any(dim(x) != dim(alpha))) 
      stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
    pd <- vector(length = nrow(x))
    for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
                                                           ])
    pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
    pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
    return(pd)
}

#decay.spec
decay.spec <- function (qm, p, lambda) 
{
  ld <- seq(1:p)^-lambda
  if (qm == 12) {
    j <- ceiling(p/3)^-lambda
    b <- 0
    if (p > 1) {
      b <- (log(1) - log(j))/(1 - p)
    }
    a <- exp(-b)
    ld <- a * exp(b * seq(1:p))
  }
  (ts(ld))
}

#dfev
dfev <- function (varobj, A0 = NULL, k) 
{
  if (!(k > 0)) {
    stop("argument 'k' in dfev() must be greater than 0")
  }
  else if (class(varobj) == c("VAR") || class(varobj) == c("BVAR")) {
    return(dfev.VAR(varobj, A0 = t(chol(varobj$mean.S)), 
                    k))
  }
  else if (class(varobj) == c("BSVAR")) {
    return(dfev.VAR(varobj, A0 = solve(varobj$A0.mode), k))
  }
}

#forc.ecdf
forc.ecdf <- function (forecasts, probs = c(0.05, 0.95), start = c(0, 1), 
          ...) 
{
  m.forecast <- apply(forecasts, 2, mean)
  quant <- apply(forecasts, 2, quantile, probs)
  vplus <- quant[1, ]
  vminus <- quant[2, ]
  forc.ci <- ts(t(rbind(m.forecast, vplus, vminus)), start = start, 
                ...)
  attr(forc.ci, "class") <- c("forc.ecdf", "mts", "ts")
  return(forc.ci)
}

#forecast
forecast <- function (varobj, nsteps, A0 = t(chol(varobj$mean.S)), shocks = matrix(0, 
                                                                                   nrow = nsteps, ncol = dim(varobj$ar.coefs)[1]), exog.fut = matrix(0, 
                                                                                                                                                     nrow = nsteps, ncol = nrow(varobj$exog.coefs)), N1, N2) 
{
  if (inherits(varobj, "VAR")) {
    return(forecast.VAR(varobj, nsteps, A0 = A0, shocks = shocks, 
                        exog.fut = exog.fut))
  }
  if (inherits(varobj, "BVAR")) {
    return(forecast.VAR(varobj, nsteps, A0 = A0, shocks = shocks, 
                        exog.fut = exog.fut))
  }
  if (inherits(varobj, "BSVAR")) {
    return(forecast.VAR(varobj, nsteps, A0 = solve(varobj$A0.mode), 
                        shocks = shocks, exog.fut = exog.fut))
  }
  if (inherits(varobj, "MSBVAR")) {
    return(forecast.MSBVAR(x = varobj, k = nsteps, N1, N2))
  }
}

#gibbs.A0
gibbs.A0 <- function (varobj, N1, N2, thin = 1, normalization = "DistanceMLA") 
{
  tmp <- sanity.check.gibbs(list(N1 = N1, N2 = N2, thin = thin, 
                                 normalization = normalization))
  methodlist <- c("DistanceMLA", "DistanceMLAhat", "Euclidean", 
                  "PositiveDiagA", "PositiveDiagAinv")
  if (tmp) {
    method <- which(methodlist == normalization) - 1
  }
  else {
    method <- which(methodlist == tmp) - 1
  }
  cat("Normalization Method: ", normalization, "(", method, 
      ")\n")
  tmp2 <- .Call("gibbsA0.cpp", varobj, as.integer(N1), as.integer(N2), 
                as.integer(thin), as.integer(method), gibbs.setup.bsvar(varobj)$UT)
  gc()
  gc()
  class(tmp2) <- c("gibbs.A0")
  return(tmp2)
}

#gibbs.msbvar
gibbs.msbvar <- function (x, N1 = 1000, N2 = 1000, permute = TRUE, Beta.idx = NULL, 
          Sigma.idx = NULL, Q.method = "MH") 
{
  if (permute == TRUE & is.null(Beta.idx) == FALSE) {
    cat("You have set permute=TRUE which means you requested a random permutation sample.\n Beta.idx argument will be ignored")
  }
  if (permute == TRUE & is.null(Sigma.idx) == FALSE) {
    cat("You have set permute=TRUE which means you requested a random permutation sample.\n Sigma.idx argument will be ignored")
  }
  if (permute == FALSE & (is.null(Beta.idx) == TRUE & is.null(Sigma.idx) == 
                          TRUE)) {
    cat("You have set permute=FALSE, but failed to provide an identification to label the regimes.\n Set Beta.idx or Sigma.idx != NULL")
  }
  if (Q.method == "Gibbs") 
    Qsampler <- Q.drawGibbs
  if (Q.method == "MH") 
    Qsampler <- Q.drawMH
  init.model <- x$init.model
  hreg <- x$hreg
  Q <- x$Q
  fp <- x$fp
  m <- x$m
  p <- x$p
  h <- x$h
  alpha.prior <- x$alpha.prior
  TT <- nrow(fp)
  e <- hreg$e
  Sigmai <- hreg$Sigma
  Betai <- hreg$Bk
  for (j in 1:N1) {
    ss <- SS.ffbs(e, TT, m, p, h, Sigmai, Q)
    Q <- Qsampler(Q, ss$transitions, prior = alpha.prior, 
                  h)
    hreg <- hregime.reg2(h, m, p, ss$SS, init.model)
    Sout <- Sigma.draw(m, h, ss, hreg, Sigmai)
    Sigmai <- Sout$Sigmai
    Bout <- Beta.draw(m, p, h, Sigmai, hreg, init.model, 
                      Betai)
    Betai <- Bout$Betai
    e <- residual.update(m, h, init.model, Betai, e)
    pf.obj <- PermuteFlip(x = list(Betai = Betai, Sigmai = Sigmai, 
                                   ss = ss, e = e, Q = Q, df = hreg$df), h, permute, 
                          Beta.idx, Sigma.idx)
    Betai <- pf.obj$Betai
    Sigmai <- pf.obj$Sigmai
    ss <- pf.obj$ss
    e <- pf.obj$e
    Q <- pf.obj$Q
    if (j%%1000 == 0) 
      cat("Burn-in iteration : ", j, "\n")
  }
  gc()
  gc()
  ss.storage <- vector("list", N2)
  transition.storage <- array(NA, c(h, h, N2))
  Beta.storage <- matrix(NA, N2, (m^2 * p + m) * h)
  Sigma.storage <- matrix(NA, N2, m * (m + 1) * 0.5 * h)
  Q.storage <- matrix(NA, N2, h^2)
  llf <- matrix(NA, N2, 1)
  if (permute == TRUE) {
    Beta.cpm <- matrix(NA, N2, (m^2 * p + m) * h)
    Sigma.wishscale <- matrix(NA, N2, 0.5 * m * (m + 1) * 
                                h)
    tmpdim <- m * ((m * p) + 1)
    Beta.cpprec <- Beta.cpv <- matrix(NA, N2, 0.5 * tmpdim * 
                                        (tmpdim + 1) * h)
    df.storage <- matrix(NA, N2, h)
    Q.cp <- matrix(NA, N2, h^2)
  }
  for (j in 1:N2) {
    ss <- SS.ffbs(e, TT, m, p, h, Sigmai, Q)
    Q <- Qsampler(Q, ss$transitions, prior = alpha.prior, 
                  h)
    hreg <- hregime.reg2(h, m, p, ss$SS, init.model)
    Sout <- Sigma.draw(m, h, ss, hreg, Sigmai)
    Sigmai <- Sout$Sigmai
    Bout <- Beta.draw(m, p, h, Sigmai, hreg, init.model, 
                      Betai)
    Betai <- Bout$Betai
    e <- residual.update(m, h, init.model, Betai, e)
    pf.obj <- PermuteFlip(x = list(Betai = Betai, Sigmai = Sigmai, 
                                   ss = ss, e = e, Q = Q, df = hreg$df), h, permute, 
                          Beta.idx, Sigma.idx)
    Betai <- pf.obj$Betai
    Sigmai <- pf.obj$Sigmai
    ss <- pf.obj$ss
    e <- pf.obj$e
    Q <- pf.obj$Q
    df <- pf.obj$df
    Sigma.storage[j, ] <- as.vector(apply(Sigmai, 3, vech))
    Beta.storage[j, ] <- as.vector(Betai)
    ss.storage[[j]] <- as.bit.integer(as.integer(ss$SS[, 
                                                       1:(h - 1)]))
    transition.storage[, , j] <- ss$transitions
    Q.storage[j, ] <- as.vector(Q)
    llf[j, 1] <- ss$llf
    if (permute == TRUE) {
      Beta.cpm[j, ] <- as.vector(Bout$Beta.cpm)
      Beta.cpprec[j, ] <- as.vector(apply(Bout$Beta.cpprec, 
                                          3, vech))
      Beta.cpv[j, ] <- as.vector(apply(Bout$Beta.cpv, 3, 
                                       vech))
      Sigma.wishscale[j, ] <- as.vector(apply(Sout$wishscale, 
                                              3, vech))
      df.storage[j, ] <- df
      Q.cp[j, ] <- as.vector(ss$transitions + alpha.prior)
    }
    if (j%%1000 == 0) 
      cat("Final iteration : ", j, "\n")
  }
  gc()
  gc()
  class(ss.storage) <- c("SS")
  if (permute == FALSE) {
    output <- list(Beta.sample = mcmc(Beta.storage), Sigma.sample = mcmc(Sigma.storage), 
                   Q.sample = mcmc(Q.storage), transition.sample = transition.storage, 
                   ss.sample = ss.storage, llf = llf, init.model = init.model, 
                   alpha.prior = alpha.prior, h = h, p = p, m = m)
  }
  else {
    output <- list(Beta.sample = mcmc(Beta.storage), Sigma.sample = mcmc(Sigma.storage), 
                   Q.sample = mcmc(Q.storage), transition.sample = transition.storage, 
                   ss.sample = ss.storage, llf = llf, init.model = init.model, 
                   alpha.prior = alpha.prior, h = h, p = p, m = m, Beta.cpm = Beta.cpm, 
                   Beta.cpprec = Beta.cpprec, Beta.cpv = Beta.cpv, Sigma.wishscale = Sigma.wishscale, 
                   df = df.storage, Q.cp = Q.cp)
  }
  class(output) <- c("MSBVAR")
  attr(output, "eqnames") <- attr(init.model, "eqnames")
  attr(output, "permute") <- permute
  attr(output, "Beta.idx") <- Beta.idx
  attr(output, "Sigma.idx") <- Sigma.idx
  attr(output, "Qsampler") <- Q.method
  tmp <- tsp(init.model$y)
  attr(output, "start") <- tmp[1]
  attr(output, "end") <- tmp[2]
  attr(output, "freq") <- tmp[3]
  return(output)
}

#granger.test
granger.test <- function (y, p) 
{
  m <- ncol(y)
  if (m < 2) {
    stop(paste("Error: granger.test needs at least 2 variables"))
  }
  results <- matrix(0, m * (m - 1), 2)
  namelist <- vector(mode = "character", m * (m - 1))
  varnames <- dimnames(y)[[2]]
  k <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      if (i == j) {
        next
      }
      Y <- embed(cbind(y[, i], y[, j]), p + 1)
      X1 <- Y[, -(1:2)]
      X2 <- X1[, ((1:p) * 2) - (1%%2)]
      restricted <- lm(Y[, 1] ~ X2)
      unrestricted <- lm(Y[, 1] ~ X1)
      ssqR <- sum(restricted$resid^2)
      ssqU <- sum(unrestricted$resid^2)
      ftest <- ((ssqR - ssqU)/p)/(ssqU/(nrow(Y) - 2 * p - 
                                          1))
      k <- k + 1
      endog.name <- varnames[i]
      exog.name <- varnames[j]
      name <- paste(exog.name, "->", endog.name)
      namelist[k] <- name
      results[k, ] <- c(ftest, 1 - pf(ftest, p, nrow(Y) - 
                                        2 * p - 1))
    }
  }
  rownames(results) <- namelist
  colnames(results) <- c("F-statistic", "p-value")
  return(results)
}

#hc.forecast
hc.forecast <- function (varobj, yconst, nsteps, burnin, gibbs, exog = NULL) 
{
  if (inherits(varobj, "VAR")) {
    stop("Not yet implemented for VAR models!\nUse a BVAR model.")
  }
  if (inherits(varobj, "BVAR")) {
    output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin, 
                              gibbs, exog)
    attr(output, "class") <- c("forecast.VAR")
    return(output)
  }
  if (inherits(varobj, "BSVAR")) {
    stop("Not yet implemented for B-SVAR models!\n")
  }
}

#initialize.msbvar
initialize.msbvar <- function (y, p, z = NULL, lambda0, lambda1, lambda3, lambda4, 
          lambda5, mu5, mu6, nu = NULL, qm, prior, h, Q = NULL) 
{
  m <- ncol(y)
  if (is.null(nu)) 
    nu <- m
  tmp <- szbvar(y, p, z = z, lambda0, lambda1, lambda3, lambda4, 
                lambda5, mu5, mu6, nu = m, qm = qm, prior = prior, posterior.fit = FALSE)
  thetahat.start <- array(NA, c(m, 1 + m * p + m, h))
  regimes <- kmeans(tmp$residuals[(m + 1):nrow(tmp$residuals), 
                                  ], centers = h)$cluster
  fptmp <- regimes[(p + 1):length(regimes)]
  df <- table(fptmp)
  for (i in 1:h) {
    s <- c(rep(i, m + 1), fptmp)
    Y1 <- tmp$Y[s == i, ]
    X1 <- tmp$X[s == i, ]
    XX <- crossprod(X1) + tmp$H0
    reg <- qr.coef(qr(crossprod(X1)), crossprod(X1, Y1) + 
                     tmp$H0[, 1:m])
    e <- Y1 - X1 %*% reg
    cSigma <- (tmp$S0 + tmp$H0[1:m, 1:m] + crossprod(e))/(df[i] + 
                                                            nu)
    thetahat.start[1:m, 1, i] <- (reg[((m * p) + 1), ])
    thetahat.start[1:m, 2:((m * p) + 1), i] <- t(reg[1:(m * 
                                                          p), ])
    thetahat.start[1:m, (2 + (m * p)):ncol(thetahat.start), 
                   i] <- t(cSigma)
  }
  if (is.null(Q) == TRUE) {
    lt <- length(fptmp)
    Qtmp <- table(fptmp[2:lt], fptmp[1:(lt - 1)])
    Q <- Qtmp/rowSums(Qtmp)
  }
  if (is.null(Q) == FALSE) {
    if (sum(ifelse(rowSums(Q) == 1, 1, 0)) < h) {
      stop("initialize.msbvar(): Invalid user input: Improper initial Q, transition matrix, for the MS process.  Rows must sum to 1.\n")
    }
  }
  return(list(init.model = tmp, thetahat.start = thetahat.start, 
              Qhat.start = Q))
}

#irf
irf <- function (varobj, nsteps, A0 = NULL) {
  if (inherits(varobj, "VAR")) {
    return(irf.VAR(varobj, nsteps, A0 = chol(varobj$mean.S)))
  }
  if (inherits(varobj, "BVAR")) {
    return(irf.BVAR(varobj, nsteps, A0 = chol(varobj$mean.S)))
  }
  if (inherits(varobj, "BSVAR")) {
    return(irf.BSVAR(varobj, nsteps, A0 = solve(varobj$A0.mode)))
  }
}

#ldwishart
ldwishart <- function (W, v, S) 
{
  k <- nrow(S)
  lgammapart <- 0
  for (i in 1:k) {
    lgammapart <- lgammapart + lgamma((v + 1 - i)/2)
  }
  denom <- lgammapart + (v * k/2) * log(2) + (k * (k - 1)/4) * 
    log(pi)
  detS <- determinant(S)$modulus[1]
  detW <- determinant(W)$modulus[1]
  hold <- solve(S) %*% W
  tracehold <- sum(hold[row(hold) == col(hold)])
  num <- detS * (-v/2) + detW * ((v - k - 1)/2) + (-1/2 * tracehold)
  return(num - denom)
}

#list.print
list.print <- function (x) 
{
  if (is.list(x$values)) {
    cat("==========================================\n")
    cat(x$labels[1], "\n")
    cat("==========================================\n")
    for (i in 1:length(x$values)) list.print(x$values[[i]])
  }
  else {
    if (length(dim(x$values)) == 3) {
      cat(x$labels[1], ": \n", sep = "")
      for (j in 1:dim(x$values)[3]) {
        cat("B(", j, ")\n", sep = "")
        prmatrix(x$values[, , j])
        cat("\n")
      }
    }
    else if (length(dim(x$values)) == 2) {
      cat(x$labels[1], ": \n", sep = "")
      prmatrix(x$values)
    }
    else if (is.null(dim(x$values))) {
      if ((length(x$values)/length(x$labels)) == 1) {
        for (i in 1:length(x$values)) cat(x$labels[i], 
                                          ":    ", x$values[i], "\n")
      }
      else {
        cat(x$labels[1], ":\n")
        for (i in 1:length(x$values)) cat(x$values[i], 
                                          "\t")
        cat("\n")
      }
    }
    cat("------------------------------------------\n")
  }
}

#mae
mae <- function (m1, m2) 
{
  tmp <- mean(abs(m1 - m2))
  return(tmp)
}

#mc.irf
mc.irf <- function (varobj, nsteps, draws = 1000, A0.posterior = NULL, 
          sign.list = rep(1, ncol(varobj$Y))) 
{
  if (inherits(varobj, "VAR")) {
    return(mc.irf.VAR(varobj = varobj, nsteps = nsteps, draws = draws))
  }
  if (inherits(varobj, "BVAR")) {
    return(mc.irf.BVAR(varobj = varobj, nsteps = nsteps, 
                       draws = draws))
  }
  if (inherits(varobj, "BSVAR")) {
    return(mc.irf.BSVAR(varobj = varobj, nsteps = nsteps, 
                        A0.posterior = A0.posterior, sign.list = sign.list))
  }
  if (inherits(varobj, "MSBVAR")) {
    return(mc.irf.MSBVAR(varobj = varobj, nsteps = nsteps, 
                         draws = length(varobj$ss.sample)))
  }
}

#mcmc.szbsvar
mcmc.szbsvar <- function (varobj, A0.posterior) 
{
  m <- dim(varobj$ar.coefs)[1]
  p <- dim(varobj$ar.coefs)[3]
  ncoef <- dim(varobj$B.posterior)[1]
  n0 <- varobj$n0
  n0cum <- c(0, cumsum(n0))
  N2 <- A0.posterior$N2
  XXinv <- chol(solve(varobj$Hpinv.posterior[[1]]))
  B.sample <- matrix(0, nrow = N2, ncol = ncoef * m)
  for (i in 1:N2) {
    A0 <- A0.get(A0.posterior$A0.posterior, i)
    A0inv <- solve(A0)
    bj <- a2b(A0, varobj$Ui)
    F.draw <- matrix(0, ncoef, m)
    for (j in 1:m) {
      btmp <- bj[(n0cum[j] + 1):(n0cum[(j + 1)])]
      F.draw[, j] <- varobj$P.posterior[[j]] %*% (btmp)
    }
    F.draw <- F.draw + XXinv %*% matrix(rnorm(m * ncoef), 
                                        ncoef, m)
    B.draw <- F.draw %*% (A0inv)
    B.sample[i, ] <- matrix(B.draw, 1, ncoef * m)
    if (i%%1000 == 0) {
      cat("Monte Carlo Iteration = ", i, "\n")
    }
  }
  output <- list(B.sample = B.sample)
  class(output) <- c("mcmc.posterior.BSVAR")
  return(output)
}

#mean.SS
mean.SS <- function (x, ...) 
{
  sums <- sum.SS(x)
  N <- length(x$ss.sample)
  return(sums/N)
}

#mountains
mountains <- function (fcasts1, fcasts2, varnames, pts, ...) 
{
  bw.1 <- dpill(fcasts1[, 1], fcasts1[, 2])
  bw.2 <- dpill(fcasts2[, 1], fcasts2[, 2])
  hill1 <- bkde2D(fcasts1, bandwidth = bw.1)
  hill2 <- bkde2D(fcasts2, bandwidth = bw.2)
  rangex <- range(hill1$x1, hill2$x1)
  rangey <- range(hill1$x2, hill2$x2)
  slice1.1 <- bkde(fcasts1[, 1], bandwidth = bw.1)
  slice1.2 <- bkde(fcasts1[, 2], bandwidth = bw.1)
  slice2.1 <- bkde(fcasts2[, 1], bandwidth = bw.2)
  slice2.2 <- bkde(fcasts2[, 2], bandwidth = bw.2)
  par(mai = c(0.75, 0.75, 0.2, 0.2), cex = 0.5)
  nf <- layout(matrix(c(1, 3, 2, 4), 2, 2, byrow = T))
  plot(slice1.1, type = "l", col = c("black"), xlim = range(slice1.1$x, 
                                                            slice2.1$x, pts[1]), ylim = range(slice1.1$y, slice2.1$y), 
       xlab = varnames[1], ylab = "Density", lwd = 2)
  lines(slice2.1, col = c("red"), type = "l", lwd = 2)
  abline(v = pts[1])
  plot(slice1.2, type = "l", col = c("black"), xlim = range(slice1.2$x, 
                                                            slice2.2$x, pts[2]), ylim = range(slice1.2$y, slice2.2$y), 
       xlab = varnames[2], ylab = "Density", lwd = 2)
  lines(slice2.2, col = c("red"), type = "l", lwd = 2)
  abline(v = pts[2])
  contour(hill1$x1, hill1$x2, hill1$fhat, col = c("black"), 
          xlim = range(rangex, pts[1]), ylim = range(rangey, pts[2]), 
          drawlabels = F)
  par(new = T)
  contour(hill2$x1, hill2$x2, hill2$fhat, col = c("red"), xlim = range(rangex, 
                                                                       pts[1]), ylim = range(rangey, pts[2]), xlab = varnames[1], 
          ylab = varnames[2], drawlabels = F)
  points(pts[1], pts[2], pch = 19)
  points(pts[3], pts[4], pch = 8)
  par(mai = c(0.5, 0.25, 0.1, 0.25), cex = 0.4)
  persp(hill1$x1, hill1$x2, hill1$fhat, xlim = rangex, ylim = rangey, 
        zlim = range(hill1$fhat, hill2$fhat), shade = 0.3, xlab = "", 
        ylab = "", zlab = "", ...)
  par(new = TRUE)
  persp(hill2$x1, hill2$x2, hill2$fhat, xlim = rangex, ylim = rangey, 
        zlim = range(hill1$fhat, hill2$fhat), ticktype = "detailed", 
        border = c("red"), shade = 0.3, cex.axis = 0.5, xlab = varnames[1], 
        ylab = varnames[2], zlab = "Density", ...)
}

#msbvar
msbvar <- function (Y, z = NULL, p, h, lambda0, lambda1, lambda3, lambda4, 
          lambda5, mu5, mu6, qm, alpha.prior = 100 * diag(h) + matrix(2, 
                                                                      h, h), prior = 0, max.iter = 40, initialize.opt = NULL) 
{
  m <- ncol(Y)
  if (h == 1) {
    stop("\n\n\t -- For MSBVAR models, h>1.  Otherwise, just for a BVAR or VAR!\n")
  }
  chk <- dim(alpha.prior)
  if (chk[1] != h) 
    stop("Incorrect number of rows in alpha.prior.")
  if (chk[2] != h) 
    stop("Incorrect number of columns in alpha.prior.")
  if (is.null(initialize.opt) == TRUE) {
    setup <- initialize.msbvar(Y, p, z, lambda0, lambda1, 
                               lambda3, lambda4, lambda5, mu5, mu6, nu = m, qm, 
                               prior, h, Q = NULL)
    init.model <- setup$init.model
    Qhat.start <- setup$Qhat.start
    thetahat.start <- setup$thetahat.start
  }
  else {
    if (inherits(initialize.opt$init.model, "BVAR") == FALSE) {
      stop("msbvar() initialize.opt list must have an object named init.model of class BVAR.  Create this using szbvar()\n")
    }
    tmp <- dim(initialize.opt$thetahat.start)
    if (tmp[1] != m) {
      stop("initialize.opt$thetahat.start has the wrong number of rows\n")
    }
    if (tmp[2] != (1 + m * p + m)) {
      stop("initialize.opt$thetahat.start has the wrong number of columns\n")
    }
    if (tmp[3] != h) {
      stop("initialize.opt$thetahat.start has the wrong array dimension\n")
    }
    if (sum(initialize.opt$Qhat.start) != h) {
      stop("msbvar(): Improper initial Q, transition matrix, for the MS process.  Rows must sum to 1.\n")
    }
    init.model <- initialize.opt$init.model
    Qhat.start <- initialize.opt$Qhat.start
    thetahat.start <- initialize.opt$thetahat.start
    tmp <- dim(alpha.prior)
    if (tmp[1] != h) {
      stop("msbvar(): Incorrect number of rows in alpha.prior")
    }
    if (tmp[2] != h) {
      stop("msbvar(): Incorrect number of columns in alpha.prior")
    }
  }
  indms <- "IAH"
  mlemod <- blkopt(Y, p, thetahat.start, Qhat.start, niter = max.iter, 
                   indms)
  hreg <- setupGibbs(mlemod, Y, init.model)
  output <- list(init.model = init.model, hreg = hreg, Q = mlemod$Qhat, 
                 fp = mlemod$fpH, m = m, p = p, h = h, alpha.prior = alpha.prior)
  class(output) <- c("MSVARsetup")
  attr(output, "eqnames") <- colnames(Y)
  return(output)
}

#msvar
msvar <- function (Y, p, h, niterblkopt = 10) 
{
  indms <- "IAH"
  n <- nrow(Y)
  m <- ncol(Y)
  if (h < 2) 
    stop("h should be an integer >=2")
  init.model <- szbvar(ts(Y), p, lambda0 = 1, lambda1 = 1, 
                       lambda3 = 1, lambda4 = 1, lambda5 = 1, mu5 = 0, mu6 = 0, 
                       prior = 2)
  Qhat.start <- (1 - (h * 0.1/(h - 1))) * diag(h) + matrix(0.1/(h - 
                                                                  1), h, h)
  thetahat.start <- array(NA, c(m, 1 + m * p + m, h))
  thetahat.start[, 1:(1 + m * p), ] <- 0
  res.im <- init.model$residuals[(m + 2):n, ]
  sig2.start <- (1/n) * crossprod(res.im, res.im)
  for (i in 1:h) {
    thetahat.start[, (1 + m * p + 1):(1 + m * p + m), i] <- sig2.start
  }
  blkopt.est <- blkopt(Y = Y, p = p, thetahat.start = thetahat.start, 
                       Qhat.start = Qhat.start, niter = niterblkopt, indms)
  hreg <- hregime.reg2.mle(h, m, p, TT = (n - p), fp = blkopt.est$fpH, 
                           init.model)
  output <- list(init.model = init.model, hreg = hreg, Q = blkopt.est$Qhat, 
                 fp = blkopt.est$fpH, m = m, p = p, h = h, llfval = blkopt.est$llfval, 
                 DirectBFGSLastSuccess = blkopt.est$DirectBFGSLastSuccess)
  class(output) <- "MSVAR"
  return(output)
}

#normalize.svar
normalize.svar <- function (A0unnormalized, A0mode, method = c("DistanceMLA", "DistanceMLAhat", 
                                             "Euclidean", "PositiveDiagA", "PositiveDiagAinv", "Unnormalized"), 
          switch.count = 0) 
{
  A0normalized <- A0unnormalized
  if (method == "DistanceMLA") {
    indx <- which(diag(solve(A0unnormalized) %*% A0mode) < 
                    0)
    if (length(indx) != 0) {
      A0normalized[, indx] <- -A0unnormalized[, indx]
      switch.count <- switch.count + 1
    }
  }
  else if (method == "DistanceMLAhat") {
    indx <- which(diag(solve(A0mode) %*% A0unnormalized) < 
                    0)
    if (length(indx) != 0) {
      A0normalized[, indx] <- -A0unnormalized[, indx]
      switch.count <- switch.count + 1
    }
  }
  else if (method == "Euclidean") {
    Adiff <- (A0unnormalized - A0mode)^2
    Adiffn <- (-A0unnormalized - A0mode)^2
    cAdiff <- colSums(Adiff)
    cAdiffn <- colSums(Adiffn)
    indx <- which(cAdiffn < cAdiff)
    if (length(indx) != 0) {
      A0normalized[, indx] <- -A0unnormalized[, indx]
      switch.count <- switch.count + 1
    }
  }
  else if (method == "PositiveDiagA") {
    indx <- which(diag(A0unnormalized) < 0)
    if (length(indx) != 0) {
      A0normalized[, indx] <- -A0unnormalized[, indx]
      switch.count <- switch.count + 1
    }
  }
  else if (method == "PositiveDiagAinv") {
    indx <- which(diag(solve(A0unnormalized)) < 0)
    if (length(indx) != 0) {
      A0normalized[, indx] <- -A0unnormalized[, indx]
      switch.count <- switch.count + 1
    }
  }
  else if (method == "Unnormalized") {
    switch.count <- switch.count
  }
  else stop("No valid normalization rule selected.  Check 'method' argument.")
  return(list(A0normalized = A0normalized, switch.count = switch.count))
}


# null.space
null.space <- function (x) 
{
  tmp <- svd(x)
  return(tmp$v[, which(tmp$d != 0)])
}


#plot.forc.ecdf
plot.forc.ecdf <- function (x, probs = c(0.05, 0.95), xlab = "", ylab = "", ylim = NA, 
          ...) 
{
  forecasts <- x
  m.forecast <- apply(forecasts, 2, mean)
  quant <- apply(forecasts, 2, quantile, probs = probs)
  vplus <- quant[1, ] + 1
  vminus <- quant[2, ] - 1
  forc.ci <- as.ts(t(rbind(m.forecast, vplus, vminus)))
  par(las = 1)
  if (is.na(ylim) == T) {
    ts.plot(forc.ci, gpars = list(xlab = xlab, ylab = ylab, 
                                  lty = c(1, 2, 2), axes = T, ...))
  }
  else {
    ts.plot(forc.ci, gpars = list(xlab = xlab, ylab = ylab, 
                                  lty = c(1, 2, 2), axes = T, ylim = ylim, ...))
  }
  box()
}


#plot.forecast
plot.forecast <- function (x, ...) 
{
  if (inherits(x, "forecast.VAR")) {
    return(plot.forecast.VAR(x = x, ...))
  }
  if (inherits(x, "forecast.BVAR")) {
    return(plot.forecast.BVAR(x = x, ...))
  }
  if (inherits(x, "forecast.BSVAR")) {
    return(plot.forecast.BSVAR(x = x, ...))
  }
}


#plot.gibbs.A0
plot.gibbs.A0 <- function (x, hpd = 0.68, varnames = attr(x, "eqnames"), ...) 
{
  m <- ncol(x$ident)
  ident <- t(x$ident)
  prob <- hpd * 100
  x <- A02mcmc(x)
  k <- 1
  par(mar = c(2, 2, 1, 1))
  split.screen(c(m, m))
  for (i in 1:m) {
    for (j in 1:m) {
      if (ident[i, j] == 1) {
        screen((i - 1) * m + j)
        den1 <- density(x[, k])
        hdr1 <- hdr(x[, k], prob = prob, den = den1)
        if (i == 1 | j == 1) 
          par(omi = c(0.15, 0.5, 0.5, 0.15))
        plot(den1, lty = 1, main = "", ylab = "", cex.axis = 0.75)
        abline(v = 0)
        nregions <- nrow(hdr1$hdr)
        maxden1 <- max(den1$y)
        for (l in 1:nregions) {
          lines(range(den1$x), rep(hdr1$falpha[l], 2), 
                lty = 2)
          for (n in 1:length(hdr1$hdr[l, ])) lines(rep(hdr1$hdr[l, 
                                                                n], 2), c((0.01 + (l - 1) * 0.02) * maxden1, 
                                                                          hdr1$falpha[l]), lty = 2)
        }
        for (l in 1:nrow(hdr1$hdr)) add.hdr(hdr1$hdr[l, 
                                                     ], (0.01 + (l - 1) * 0.02) * maxden1, 0.1 * 
                                              maxden1, col = "black", horiz = TRUE, border = FALSE)
        if (i == 1) 
          mtext(varnames[j], side = 3, line = 1, outer = FALSE)
        if (j == 1) 
          mtext(varnames[i], side = 2, line = 3, outer = FALSE)
        k <- k + 1
      }
      if (ident[i, j] == 0 & (i == 1 || j == 1)) {
        screen((i - 1) * m + j)
        plot.new()
        if (i == 1) 
          mtext(varnames[j], side = 3, line = 1, outer = FALSE)
        if (j == 1) 
          mtext(varnames[i], side = 2, line = 3, outer = FALSE)
      }
    }
  }
  close.screen(all.screens = TRUE)
}


#plot.irf
plot.irf <- function (x, varnames = attr(x, "eqnames"), ...) 
{
  if (inherits(x, "irf.VAR")) {
    plot.irf.VAR(x, varnames = varnames, ...)
  }
  if (inherits(x, "irf.BVAR")) {
    plot.irf.BVAR(x, varnames = varnames, ...)
  }
  if (inherits(x, "irf.BSVAR")) {
    plot.irf.BSVAR(x, varnames = varnames, ...)
  }
}


#plot.mc.irf
plot.mc.irf <- function (x, method = c("Sims-Zha2"), component = 1, probs = c(0.16, 
                                                               0.84), varnames = attr(x, "eqnames"), regimelabels = NULL, 
          ask = TRUE, ...) 
{
  if (inherits(x, "mc.irf.VAR")) {
    tmp <- plot.mc.irf.VAR(x, method = method, component = component, 
                           probs = probs, varnames = varnames, ...)
  }
  if (inherits(x, "mc.irf.BVAR")) {
    tmp <- plot.mc.irf.BVAR(x, method = method, component = component, 
                            probs = probs, varnames = varnames, ...)
  }
  if (inherits(x, "mc.irf.BSVAR")) {
    tmp <- plot.mc.irf.BSVAR(x, method = method, component = component, 
                             probs = probs, varnames = varnames, ...)
  }
  if (inherits(x, "mc.irf.MSBVAR")) {
    tmp <- plot.mc.irf.MSBVAR(x, method = method, component = component, 
                              probs = probs, varnames = varnames, regimelabels = regimelabels, 
                              ask, ...)
  }
  return(invisible(tmp))
}


#plot.ms.irf
plot.ms.irf <- function (x, method = "Sims-Zha2", component = 1, probs = c(0.16, 
                                                            0.84), varnames = attr(x, "eqnames"), ...) 
{
  dd <- dim(x$shortrun)
  N2 <- dd[1]
  nsteps <- dd[2]
  h <- dd[4]
  m <- sqrt(dd[3])
  out <- vector(mode = "list", length = (h + 1))
  rgs <- array(0, c(2, m^2, h))
  for (i in 1:h) {
    out[[i]] <- compute.plot.mc.irf(x$shortrun[, , , i], 
                                    method, component, probs)
    rgs[, , i] <- apply(out[[i]]$responses, 3, range)
  }
  minmax <- matrix(0, nrow = m, ncol = 2)
  tmp <- (c(1:m^2)%%m)
  tmp[tmp == 0] <- m
  indices <- sort(tmp, index.return = T)$ix
  dim(indices) <- c(m, m)
  for (i in 1:m) {
    minmax[i, ] <- range(rgs[, indices[, i], 1:h])
  }
  col <- rep(1:h, each = 3)
  j <- 1
  par(mfcol = c(m, m), mai = c(0.25, 0.25, 0.15, 0.25), omi = c(0.15, 
                                                                0.75, 1, 0.15))
  for (i in 1:m^2) {
    lims <- ifelse((i - m)%%m == 0, m, (i - m)%%m)
    for (k in 1:h) {
      if (k == 1) {
        tmp <- out[[k]]$responses[, , i]
      }
      else {
        tmp <- cbind(tmp, out[[k]]$responses[, , i])
      }
    }
    ts.plot(tmp, gpars = list(xlab = "", ylab = "", ylim = minmax[lims, 
                                                                  ], col = col, ...))
    abline(h = 0)
    if (i <= m) {
      mtext(varnames[i], side = 2, line = 3)
    }
    if ((i - 1)%%m == 0) {
      mtext(varnames[j], side = 3, line = 2)
      j <- j + 1
    }
  }
  mtext("Response in", side = 2, line = 3, outer = T)
  mtext("Shock to", side = 3, line = 3, outer = T)
  return(invisible(out))
}


#plot.SS
plot.SS <- function (x, ylab = "State Probabilities", ...) 
{
  tmp <- mean.SS(x)
  shift <- x$p/attr(x, "freq")
  plot(ts(tmp, start = attr(x, "start") + shift, end = attr(x, 
                                                            "end"), frequency = attr(x, "freq")), plot.type = "single", 
       col = 1:ncol(tmp), ylim = c(0, 1), ylab = ylab, ...)
  abline(h = 0.5, lty = 2, ...)
}


#plotregimeid
plotregimeid <- function (x, type = c("all", "intercepts", "AR1", "Sigma", "Q"), 
          ask = TRUE, ...) 
{
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  devAskNewPage()
  m <- x$m
  p <- x$p
  h <- x$h
  N2 <- length(x$ss.sample)
  mpplus1 <- m * p + 1
  nc <- m * (mpplus1)
  ii <- seq(mpplus1, by = mpplus1, length = m)
  aii <- rep(1:m, m) + (rep(0:(m - 1), each = m)) * mpplus1
  eqnnames <- colnames(x$init.model$y)
  if (type == "all" || type == "intercepts") {
    Beta <- matrix(t(x$Beta.sample), byrow = TRUE, ncol = nc)
    intercepts <- as.data.frame(Beta[, ii])
    colnames(intercepts) <- eqnnames
    cl.int <- kmeans(intercepts, centers = h, nstart = 10)
    pairs(intercepts, pch = ".", col = cl.int$cluster)
    title("Intercepts pairs by regime", line = 3)
    form <- as.formula(paste("~", paste(lapply(names(intercepts), 
                                               as.name), collapse = "+")))
    devAskNewPage(ask = ask)
    print(densityplot(form, data = intercepts, outer = TRUE, 
                      groups = cl.int$cluster, xlab = NULL, default.scales = list(relation = "free"), 
                      plot.points = "rug", main = "Intercept densities by regime"))
    form <- eval(parse(text = paste(paste(lapply(names(intercepts), 
                                                 as.name), collapse = "+"), "~ idx")))
    idx <- 1:nrow(intercepts)
    devAskNewPage(ask = ask)
    print(xyplot(form, data = intercepts, groups = cl.int$cluster, 
                 type = "l", ylab = NULL, default.scales = list(relation = "free"), 
                 main = "Intercept traceplots by regime"))
  }
  if (type == "all" || type == "AR1") {
    Beta <- matrix(t(x$Beta.sample), byrow = TRUE, ncol = nc)
    ar1 <- as.data.frame(Beta[, aii])
    idx <- expand.grid(1:m, 1:m)
    tmp <- cbind(rep(eqnnames, m), idx)
    colnames(ar1) <- paste(tmp[, 1], "(", tmp[, 2], ",", 
                           tmp[, 3], ")", sep = "")
    cl.ar1 <- kmeans(ar1, centers = h, nstart = 10)
    form <- as.formula(paste("~", paste(lapply(names(ar1), 
                                               as.name), collapse = "+")))
    devAskNewPage(ask = ask)
    print(densityplot(form, data = ar1, outer = TRUE, groups = cl.ar1$cluster, 
                      xlab = NULL, default.scales = list(relation = "free"), 
                      plot.points = "rug", main = "AR(1) densities by regime"))
    form <- eval(parse(text = paste(paste(lapply(names(ar1), 
                                                 as.name), collapse = "+"), "~ idx")))
    idx <- 1:nrow(ar1)
    devAskNewPage(ask = ask)
    print(xyplot(form, data = ar1, groups = cl.ar1$cluster, 
                 type = "l", ylab = NULL, main = "AR(1) coefficients traceplots by regime"))
  }
  if (type == "all" || type == "Sigma") {
    nvar <- (m * (m + 1) * 0.5)
    Sigmaout <- matrix(t(x$Sigma.sample), byrow = TRUE, ncol = nvar)
    tmp <- m:2
    tmp <- rep(c(1, tmp), h)
    ssidx <- tmp
    for (i in 2:(m * h)) ssidx[i] <- ssidx[i - 1] + tmp[i]
    Sigmaout <- as.data.frame(Sigmaout[, ssidx[1:m]])
    colnames(Sigmaout) <- eqnnames
    cl.sigma <- kmeans(Sigmaout, centers = h, nstart = 10)
    devAskNewPage(ask = ask)
    pairs(Sigmaout, pch = ".", col = cl.sigma$cluster, ...)
    title("Variances pairs plot by regime", line = 3)
    form <- as.formula(paste("~", paste(lapply(names(Sigmaout), 
                                               as.name), collapse = "+")))
    devAskNewPage(ask = ask)
    print(densityplot(form, data = Sigmaout, outer = TRUE, 
                      groups = cl.sigma$cluster, xlab = NULL, default.scales = list(relation = "free"), 
                      plot.points = "rug", main = "Variance densities by regime"), 
          ...)
    form <- eval(parse(text = paste(paste(lapply(names(Sigmaout), 
                                                 as.name), collapse = "+"), "~ idx")))
    idx <- 1:nrow(Sigmaout)
    devAskNewPage(ask = ask)
    print(xyplot(form, data = Sigmaout, groups = cl.sigma$cluster, 
                 type = "l", ylab = NULL, default.scales = list(relation = "free"), 
                 main = "Variance traceplots by regime"), ...)
  }
  if (type == "all" || type == "Q") {
    Q <- as.data.frame(x$Q.sample)
    cl.Q <- kmeans(Q, centers = h, nstart = 10)
    idx <- cbind(rep(1:h, each = h), rep(1:h, times = h))
    Qnames <- paste("Q_", idx[, 1], idx[, 2], sep = "")
    colnames(Q) <- Qnames
    devAskNewPage(ask = ask)
    pairs(Q, pch = ".", col = cl.Q$cluster, ...)
    title("Transitions pairs plot by regime", line = 3)
    form <- as.formula(paste("~", paste(lapply(names(Q), 
                                               as.name), collapse = "+")))
    devAskNewPage(ask = ask)
    print(densityplot(form, data = Q, outer = TRUE, groups = cl.Q$cluster, 
                      xlab = NULL, default.scales = list(relation = "free"), 
                      plot.points = "rug", main = "Transition densities by regime"), 
          ...)
    form <- eval(parse(text = paste(paste(lapply(names(Q), 
                                                 as.name), collapse = "+"), "~ idx")))
    idx <- 1:nrow(Q)
    devAskNewPage(ask = ask)
    print(xyplot(form, data = Q, groups = cl.Q$cluster, type = "l", 
                 ylab = NULL, default.scales = list(relation = "free"), 
                 main = "Transition traceplots by regime"), ...)
    Beta <- matrix(t(x$Beta.sample), byrow = TRUE, ncol = nc)
    intercepts <- as.data.frame(Beta[, ii])
    colnames(intercepts) <- eqnnames
    qdiag <- diag(matrix(1:h^2, h, h))
    devAskNewPage(ask = ask)
    par(mfrow = c(2, round(m/2)), omi = c(0.5, 0.75, 0.75, 
                                          0.25))
    for (i in 1:m) {
      plot(intercepts[, i], matrix(unlist(Q[, qdiag]), 
                                   ncol = 1), pch = ".", col = cl.Q$cluster, xlab = names(intercepts)[i], 
           ylab = "Transition Probability Regimes")
    }
    title("Intercepts by transition probability regimes", 
          outer = TRUE, line = 1)
  }
  devAskNewPage(FALSE)
  invisible()
}


#posterior.fit
posterior.fit <- function (varobj, A0.posterior.obj = NULL, maxiterbs = 500) 
{
  if (inherits(varobj, "VAR")) {
    stop("posterior.fit() not implemented for VAR class objects since they do not have a proper prior.\n")
  }
  if (inherits(varobj, "BVAR")) {
    output <- posterior.fit.BVAR(varobj)
    attr(output, "class") <- c("posterior.fit.BVAR")
    return(output)
  }
  if (inherits(varobj, "BSVAR")) {
    output <- posterior.fit.BSVAR(varobj, A0.posterior.obj)
    attr(output, "class") <- c("posterior.fit.BSVAR")
    return(output)
  }
  if (inherits(varobj, "MSBVAR")) {
    output <- posterior.fit.MSBVAR(x = varobj, maxiterbs = maxiterbs)
    attr(output, "class") <- c("posterior.fit.MSBVAR")
    return(output)
  }
}


#print.dfev
print.dfev <- function (x, latex = F, file = NULL, ...) 
{
  dfev.obj <- x
  errors <- dfev.obj$errors
  names <- attr(dfev.obj, "eqnames")
  std.err <- dfev.obj$std.err
  k <- dim(errors)[1]
  m <- dim(errors)[2]
  if (latex == T) {
    for (i in 1:m) {
      tmp <- matrix(errors[, , i], nrow = k, ncol = m)
      tmp <- cbind(std.err[, i], tmp)
      colnames(tmp) <- c("Std. Error", names)
      if (i == 1) {
        if (is.null(file)) {
          print(xtable(tmp, digits = rep(1, ncol(tmp) + 
                                           1), caption = paste("Decomposition of Forecast Errors for a Shock to", 
                                                               names[i])), append = F, table.placement = "p")
        }
        else {
          print(xtable(tmp, digits = rep(1, ncol(tmp) + 
                                           1), caption = paste("Decomposition of Forecast Errors for a Shock to", 
                                                               names[i])), file = file, append = F, table.placement = "p")
        }
      }
      else {
        if (is.null(file)) {
          print(xtable(tmp, digits = rep(1, ncol(tmp) + 
                                           1), caption = paste("Decomposition of Forecast Errors for a Shock to", 
                                                               names[i])), append = T, table.placement = "p")
        }
        else {
          print(xtable(tmp, digits = rep(1, ncol(tmp) + 
                                           1), caption = paste("Decomposition of Forecast Errors for a Shock to", 
                                                               names[i])), file = file, append = T, table.placement = "p")
        }
      }
    }
  }
  else {
    for (i in 1:m) {
      cat(paste("Decomposition of Forecast Errors for a Shock to", 
                names[i], "\n"))
      cat("-------------------------------------------------------------\n")
      tmp <- matrix(errors[, , i], nrow = k, ncol = m)
      tmp <- cbind(std.err[, i], tmp)
      colnames(tmp) <- c("Std. Error", names)
      print(tmp)
      cat("-------------------------------------------------------------\n")
    }
  }
}


#print.posterior.fit
print.posterior.fit <- function (x, ...) 
{
  if (inherits(x, "posterior.fit.BVAR")) {
    print.posterior.fit.BVAR(x, ...)
  }
  if (inherits(x, "posterior.fit.BSVAR")) {
    print.posterior.fit.BSVAR(x, ...)
  }
  if (inherits(x, "posterior.fit.MSBVAR")) {
    print.posterior.fit.MSBVAR(x, ...)
  }
}


#rdirichlet
rdirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}


#reduced.form.var
reduced.form.var <- function (Y, p, z = NULL) 
{
  sanity.check.var(list(Y = Y, p = p, z = z))
  y <- Y
  dat <- as.matrix(Y)
  n <- nrow(dat)
  m <- ncol(dat)
  ncoef <- (m * p) + 1
  ndum <- m + 1
  capT <- n - p + ndum
  Ts <- n - p
  X <- matrix(0, nrow = Ts, ncol = ncoef)
  Y <- matrix(0, nrow = Ts, ncol = m)
  const <- matrix(1, nrow = Ts)
  X[, ncoef] <- const
  for (i in 1:p) {
    X[(1:Ts), (m * (i - 1) + 1):(m * i)] <- matrix(dat[(p + 
                                                          1 - i):(n - i), ], ncol = m)
  }
  X <- cbind(X[, 1:ncoef - 1], X[, ncoef])
  Y[1:Ts, ] <- matrix(dat[(p + 1):n, ], ncol = m)
  if (is.null(z) == FALSE) {
    X <- cbind(X, z[(p + 1):n, ])
  }
  XX <- crossprod(X)
  Bh <- solve(XX, crossprod(X, Y), tol = 1e-16)
  u <- (Y - X %*% (Bh))
  Sh <- crossprod(u)/Ts
  Sh1 <- crossprod(u)/capT
  intercept <- Bh[(m * p + 1), ]
  ar.coefs <- t(Bh[1:(m * p), ])
  dim(ar.coefs) <- c(m, m, p)
  ar.coefs <- aperm(ar.coefs, c(2, 1, 3))
  if (is.null(z) == FALSE) {
    exog.coefs <- Bh[(m * p + 2):nrow(Bh), ]
    num.exog <- ncol(z)
    z <- as.matrix(z)
  }
  else {
    exog.coefs <- NA
    num.exog <- 0
  }
  pfit <- list(capT = capT, ncoef = ncoef, num.exog = num.exog, 
               Sh1 = Sh1)
  output <- list(intercept = intercept, ar.coefs = ar.coefs, 
                 Bhat = Bh, vcv = Sh, exog.coefs = exog.coefs, residuals = u, 
                 mean.S = Sh, hstar = XX, X = X, Y = Y, y = y, pfit = pfit)
  class(output) <- c("VAR")
  attr(output, "eqnames") <- colnames(y)
  return(output)
}


#regimeSummary
regimeSummary <- function (x, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) 
{
  if (class(x) != "MSBVAR") {
    stop("This function can only be used with Gibbs / MCMC output from gibbs.msbvar()\n")
  }
  qij <- expand.grid(1:x$h, 1:x$h)
  colnames(x$Q.sample) <- paste("q_", qij[, 1], qij[, 2], sep = "")
  quant <- summary(x$Q.sample, quantiles = quantiles)
  cat("##############################################################\n")
  cat("Summary of MCMC draws for elements of the transition matrix Q\n")
  cat("##############################################################\n")
  print(quant)
  cat("##############################################################\n")
  cat("Full mean transition matrix\n")
  cat("##############################################################\n")
  tmp <- matrix(quant$statistics[, 1], x$h, x$h)
  colnames(tmp) <- rownames(tmp) <- paste("Regime", 1:x$h)
  print(tmp)
  N2 <- length(x$ss.sample)
  lrQ <- sapply(1:N2, function(i) {
    steady.Q(matrix(x$Q.sample[i, ], x$h, x$h))
  })
  rownames(lrQ) <- paste("Regime", 1:x$h)
  lrQquant <- summary(mcmc(t(lrQ)), quantiles = quantiles)
  cat("\n##############################################################\n")
  cat("Ergodic regime probabilities\n")
  cat("##############################################################\n")
  print(lrQquant)
  dur <- 1/(1 - lrQ)
  rownames(dur) <- paste("Regime", 1:x$h, "duration")
  durquant <- summary(mcmc(t(dur)), quantiles = quantiles)
  cat("##############################################################\n")
  cat("Regime durations and quantiles\n")
  cat("##############################################################\n")
  print(durquant)
  invisible(list(Q.summary = quant, lrQ = lrQquant, durations = durquant))
}


#restmtx
restmtx <- function (nsteps, m) 
{
  tmp <- c(1, rep(0, m - 1))
  R <- matrix(0, nrow = nsteps, ncol = (nsteps * m))
  for (i in 1:nsteps) {
    R[i, 1:(i * m)] <- matrix(rep(tmp, i), nrow = 1)
  }
  return(R)
}


#rmse
rmse <- function (m1, m2) 
{
  tmp <- sqrt(mean((m1 - m2)^2))
  return(tmp)
}


#rmultnorm
rmultnorm <- function (n, mu, vmat, tol = 1e-10) 
{
  p <- ncol(vmat)
  if (length(mu) != p) 
    stop(paste("mu vector is the wrong length:", length(mu)))
  vs <- La.svd(vmat)
  vsqrt <- t(t(vs$vt) %*% (t(vs$u) * sqrt(vs$d)))
  ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
  ans <- sweep(ans, 2, mu, "+")
  dimnames(ans) <- list(NULL, dimnames(vmat)[[2]])
  ans
}


#rwishart
rwishart <- function (N, df, Sigma) 
{
  p = nrow(Sigma)
  SqrtSigma <- t(chol(Sigma))
  tmp <- array(0, c(p, p, N))
  for (i in 1:N) {
    Z <- matrix(rnorm(df * p), nrow = p, ncol = df)
    ZS <- crossprod(Z, SqrtSigma)
    tmp[, , i] <- crossprod(ZS)
  }
  if (N == 1) {
    return(matrix(tmp, p, p))
  }
  else {
    return(tmp)
  }
}

#simulateMSAR
simulateMSAR <- function (bigt, Q, theta, st1, y1) 
{
  h <- ncol(Q)
  p <- ncol(theta) - 2
  if ((nrow(Q) != h) || (ncol(Q) != h)) 
    stop("Number of rows or columns does not equal h")
  if (sum(((Q >= 0) & (Q <= 1))) != h * h) 
    stop("Probabilities are not between zero and one")
  if (sum(rowSums(Q)) != h) 
    stop("Rows of transition matrix Q do not sum to unity).")
  if (nrow(theta) != h) 
    stop("Number of parameters in theta does not match number of states")
  if (ncol(theta) != (p + 2)) 
    stop("Number of parameters in theta does not match AR(p)")
  st.sim <- rep(NA, bigt)
  st.sim[1] <- st1
  for (i in 2:bigt) {
    p.cumsum <- cumsum(Q[st.sim[i - 1], ])
    u <- runif(1)
    st.sim[i] <- h - sum(u <= p.cumsum) + 1
  }
  emat <- NULL
  for (i in 1:h) {
    emat <- cbind(emat, rnorm(bigt, 0, sqrt(theta[i, p + 
                                                    2])))
  }
  y.sim <- rep(NA, bigt)
  y.sim[1:p] <- y1
  for (i in (p + 1):bigt) {
    s <- st.sim[i]
    ylag <- matrix(c(1, y.sim[(i - 1):(i - p)]))
    y.sim[i] <- crossprod(theta[s, 1:(p + 1)], ylag) + emat[i, 
                                                            s]
  }
  return(list(Y = y.sim, st = st.sim))
}

#simulateMSVAR
simulateMSVAR <- function (bigt, m, p, var.beta0, var.betas, e.vcv, Q, seed = 214) 
{
  h <- ncol(Q)
  st1 <- 1
  if ((nrow(Q) != h) || (ncol(Q) != h)) 
    stop("Number of rows or columns does not equal h")
  if (sum(((Q >= 0) & (Q <= 1))) != h * h) 
    stop("Probabilities are not between zero and one")
  if (sum(rowSums(Q)) != h) 
    stop("Rows of transition matrix Q do not sum to unity).")
  st.sim <- rep(NA, bigt)
  st.sim[1] <- st1
  for (i in 2:bigt) {
    p.cumsum <- cumsum(Q[st.sim[i - 1], ])
    u <- runif(1)
    st.sim[i] <- h - sum(u <= p.cumsum) + 1
  }
  e.rmvn <- array(NA, c(bigt, m, h))
  for (i in 1:h) {
    e.rmvn[, , i] <- rmvnorm(bigt, rep(0, m), e.vcv[, , i])
  }
  Y.sim <- matrix(NA, bigt, m)
  Y.sim[1:p, ] <- matrix(rep(var.beta0[, , st1], p), p, byrow = TRUE)
  for (i in p:(bigt - 1)) {
    s <- st.sim[i]
    Ylag <- matrix(c(t(Y.sim[i:(i - p + 1), ])), 1)
    Y.sim[i + 1, ] <- var.beta0[, , s] + tcrossprod(var.betas[, 
                                                              , s], Ylag) + c(e.rmvn[i, , s])
  }
  return(list(Y = Y.sim, st = st.sim))
}

#SS.ffbs
SS.ffbs <- function (e, bigt, m, p, h, sig2, Q) 
{
  ff <- .Fortran("ForwardFilter", nvar = as.integer(m), e = e, 
                 bigK = as.integer(h), bigT = as.integer(bigt), nbeta = as.integer(1 + 
                                                                                     m * p), sig2, Q, llh = double(1), pfilt = matrix(0, 
                                                                                                                                      bigt + 1, h))
  rvu <- runif(bigt + 1)
  bs <- .Fortran("BackwardSampler", bigK = as.integer(h), bigT = as.integer(bigt), 
                 ff$pfilt, Q, rvu, backsamp.bigS = integer(bigt + 1), 
                 transmat = matrix(as.integer(0), h, h))
  matSS <- matrix(as.integer(matrix(seq(1, h), bigt, h, byrow = TRUE) == 
                               bs$backsamp.bigS[-1]), bigt, h)
  ss <- list(SS = matSS, transitions = bs$transmat, llf = ff$llh)
  return(ss)
}


sum.SS <- function (x, ...) 
{
  h <- x$h
  N <- length(x$ss.sample)
  ss <- x$ss.sample
  TTh <- virtual(x$ss.sample[[1]])$Length
  TT <- TTh/(h - 1)
  sums <- apply(matrix(unlist(lapply(x$ss.sample, as.integer)), 
                       nrow = TTh, ncol = N), 1, sum)
  sums <- matrix(sums, TT, h - 1)
  sums <- cbind(sums, rep(N, TT) - rowSums(sums))
  return(sums)
}


#summary.dfev
summary.dfev <- function (object, latex = F, file = NULL, ...) 
{
  print.dfev(object, latex = F, file = NULL, ...)
}


#summary.forecast
summary.forecast <- function (object, probs = c(0.16, 0.84), ...) 
{
  if (inherits(object, "forecast.VAR")) {
    return(summary.forecast.VAR(object, probs = probs))
  }
  if (inherits(object, "forecast.BVAR")) {
    return(summary.forecast.BVAR(object, probs = probs))
  }
  if (inherits(object, "forecast.BSVAR")) {
    return(summary.forecast.BSVAR(object, probs = probs))
  }
}


#SZ.prior.evaluation
SZ.prior.evaluation <- function (Y, p, lambda0, lambda1, lambda3, lambda4, lambda5, 
          mu5, mu6, z = NULL, nu = ncol(Y) + 1, qm, prior = 0, nsteps, 
          y.future) 
{
  combos <- length(lambda0) * length(lambda1) * length(lambda3) * 
    length(lambda4) * length(lambda5) * length(mu5) * length(mu6)
  results <- matrix(0, combos, 11)
  results[, 1:7] <- as.matrix(expand.grid(lambda0 = lambda0, 
                                          lambda1 = lambda1, lambda3 = lambda3, lambda4 = lambda4, 
                                          lambda5 = lambda5, mu5 = mu5, mu6 = mu6))
  for (i in 1:nrow(results)) {
    fit <- szbvar(Y, p, z, lambda0 = results[i, 1], lambda1 = results[i, 
                                                                      2], lambda3 = results[i, 3], lambda4 = results[i, 
                                                                                                                     4], lambda5 = results[i, 5], mu5 = results[i, 6], 
                  mu6 = results[i, 7], nu = ncol(Y) + 1, qm = qm, prior = prior, 
                  posterior.fit = T)
    forecast <- forecast(fit, nsteps)
    eval.forecasts <- cf.forecasts(forecast[(nrow(forecast) - 
                                               nsteps + 1):nrow(forecast), ], y.future)
    tmp <- c(eval.forecasts[1], eval.forecasts[2], fit$marg.llf[1], 
             fit$marg.post[1])
    results[i, 8:11] <- tmp
    if (i%%100 == T) {
      cat("Finished model", i, "of", nrow(results), "\n")
    }
  }
  colnames(results) <- c("lambda0", "lambda1", "lambda3", "lambda4", 
                         "lambda5", "mu5", "mu6", "RMSE", "MAE", "LLF", "logMDD")
  return(results)
}


#szbsvar
szbsvar <- function (Y, p, z = NULL, lambda0, lambda1, lambda3, lambda4, 
          lambda5, mu5, mu6, ident, qm = 4) 
{
  sanity.check.bsvar(list(Y = Y, p = p, z = z, lambda0 = lambda0, 
                          lambda1 = lambda1, lambda3 = lambda3, lambda4 = lambda4, 
                          lambda5 = lambda5, mu5 = mu5, mu6 = mu6, qm = qm, ident = ident))
  m <- ncol(Y)
  nexog <- ifelse(is.null(z) == TRUE, 0, ncol(z))
  ncoef <- m * p + nexog + 1
  n <- nrow(Y)
  if (dim(ident)[1] != m) {
    stop("Identification matrix 'ident' and dimension of 'Y' are nonconformable.")
  }
  endog.ncoef <- (m * p) + 1
  ndum <- m + 1
  capT <- n - p + ndum
  Ts <- n - p
  Q <- array(0, c(m, m, m))
  for (i in 1:m) {
    Q[, , i] <- diag(ident[, i])
  }
  Ui <- sapply(1:m, function(i) {
    null.space(Q[, , i])
  }, simplify = F)
  Pi <- matrix(0, ncoef, m)
  diag(Pi) <- 1
  s2i <- matrix(0, nrow = m, ncol = 1)
  for (i in 1:m) {
    s2i[i, 1] <- ar.ols(Y[, i], aic = FALSE, order.max = p, 
                        intercept = TRUE, demean = FALSE)$var.pred
  }
  S0 <- diag(m)
  diag(S0) <- 1/s2i
  if (qm == 12) {
    j <- ceiling(p/3)^-lambda3
    b <- 0
    if (p > 1) {
      b <- (log(1) - log(j))/(1 - p)
    }
    a <- exp(-b)
  }
  Aplus.prior.cov <- matrix(0, ncoef, 1)
  for (i in 1:p) {
    if (qm == 12) {
      ld <- a * exp(b * i * lambda3)
    }
    for (j in 1:m) {
      if (qm == 12) {
        Aplus.prior.cov[((i - 1) * m + j), 1] <- ld^2/s2i[j, 
                                                          1]
      }
      else {
        Aplus.prior.cov[((i - 1) * m + j), 1] <- (1/i^lambda3)^2/s2i[j, 
                                                                     1]
      }
    }
  }
  A0.prior <- lambda0^2/s2i
  Aplus.prior <- lambda0^2 * lambda1^2 * Aplus.prior.cov
  Aplus.prior[(m * p + 1), 1] <- (lambda0 * lambda4)^2
  if (nexog > 0) {
    Aplus.prior[(m * p + 2):ncoef, 1] <- lambda0^2 * lambda5^2
  }
  Aplus.prior1 <- Aplus.prior[(m + 1):ncoef, 1]
  Hptd <- diag(as.vector(Aplus.prior))
  Hptdi <- diag(as.vector(1/Aplus.prior))
  H0multi <- array(0, c(m, m, m))
  H0invmulti <- H0multi
  Hpmulti <- array(0, c(ncoef, ncoef, m))
  Hpmultiinv <- Hpmulti
  H0td <- matrix(0, m, m)
  H0tdi <- H0td
  for (i in 1:m) {
    A0i <- A0.prior
    A0i.inv <- 1/A0i
    diag(H0td) <- A0i
    diag(H0tdi) <- 1/A0i
    H0multi[, , i] <- H0td
    H0invmulti[, , i] <- H0tdi
    Hpmulti[, , i] <- Hptd
    Hpmultiinv[, , i] <- Hptdi
  }
  Hpinv.tilde <- Hpmultiinv
  Pi.tilde <- sapply(1:m, function(i) {
    Pi %*% Ui[[i]]
  }, simplify = F)
  H0inv.tilde <- sapply(1:m, function(i) {
    t(Ui[[i]]) %*% H0invmulti[, , i] %*% Ui[[i]]
  }, simplify = F)
  if (is.null(z)) {
    num.exog <- 0
  }
  else {
    num.exog <- ncol(z)
    z <- as.matrix(z)
    if (det(crossprod(cbind(rep(1, nrow(z)), z))) <= 0) {
      stop("Matrix of exogenous variables, z has deficient rank.")
    }
  }
  if (p == 1) {
    datint <- as.vector(Y[1, ])
  }
  else {
    datint <- as.vector(apply(Y[1:p, ], 2, mean))
  }
  X1 <- matrix(0, nrow = capT, ncol = endog.ncoef)
  Y1 <- matrix(0, nrow = capT, ncol = m)
  const <- matrix(1, nrow = capT)
  const[1:m, ] <- 0
  X1[, endog.ncoef] <- const
  for (i in 1:m) {
    Y1[ndum, i] <- datint[i]
    Y1[i, i] <- datint[i]
    for (j in 1:p) {
      X1[ndum, m * (j - 1) + i] <- datint[i]
      X1[i, m * (j - 1) + i] <- datint[i]
    }
  }
  for (i in 1:p) {
    X1[(ndum + 1):capT, (m * (i - 1) + 1):(m * i)] <- matrix(Y[(p + 
                                                                  1 - i):(n - i), ], ncol = m)
  }
  if (is.null(z) == F) {
    pad.z <- matrix(0, nrow = capT, ncol = ncol(z))
    pad.z[(ndum + 1):capT, ] <- matrix(z[(p + 1):n, ], ncol = ncol(z))
    X1 <- cbind(X1, pad.z)
  }
  Y1[(ndum + 1):capT, ] <- matrix(Y[(p + 1):n, ], ncol = m)
  X1[1:m, ] <- mu5 * X1[1:m, ]
  Y1[1:m, ] <- mu5 * Y1[1:m, ]
  X1[ndum, ] <- mu6 * X1[ndum, ]
  Y1[ndum, ] <- mu6 * Y1[ndum, ]
  XX <- crossprod(X1)
  XY <- crossprod(X1, Y1)
  YY <- crossprod(Y1)
  Hpinv.posterior <- sapply(1:m, function(i) {
    XX + Hpinv.tilde[, , i]
  }, simplify = F)
  P1.posterior <- sapply(1:m, function(i) {
    XY %*% Ui[[i]] + Hpinv.tilde[, , i] %*% Pi.tilde[[i]]
  }, simplify = F)
  P.posterior <- sapply(1:m, function(i) {
    solve(Hpinv.posterior[[i]]) %*% P1.posterior[[i]]
  }, simplify = F)
  H0inv.posterior <- sapply(1:m, function(i) {
    (t(Ui[[i]]) %*% YY %*% Ui[[i]] + H0inv.tilde[[i]] + t(Pi.tilde[[i]]) %*% 
       Hpinv.tilde[, , i] %*% Pi.tilde[[i]] - t(P1.posterior[[i]]) %*% 
       P.posterior[[i]])
  }, simplify = F)
  n0 <- sapply(1:m, function(i) {
    ncol(as.matrix(Ui[[i]]))
  })
  n0cum <- c(0, cumsum(n0))
  b <- (1/max(s2i)) * (rnorm(sum(n0)))
  cat("Estimating starting values for the numerical optimization\nof the log posterior of A(0)\n")
  max.obj <- optim(b, A0.llf, method = c("Nelder-Mead"), control = list(maxit = 6000, 
                                                                        fnscale = capT, trace = 0), Ui = Ui, df = capT, H0inv.posterior = H0inv.posterior)
  cat("Estimating the final values for the numerical optimization\nof the log posterior of A(0)\n")
  max.obj <- optim(max.obj$par, A0.llf, method = c("BFGS"), 
                   hessian = F, control = list(maxit = 5000, fnscale = capT, 
                                               trace = 1), Ui = Ui, df = capT, H0inv.posterior = H0inv.posterior)
  if (max.obj$convergence != 0) {
    stop("Estiamtes of A(0) did not converge.  You should restart the function with a new seed.")
  }
  A0.mode <- b2a(max.obj$par, Ui)
  F.posterior <- matrix(0, ncoef, m)
  for (i in 1:m) {
    bj <- max.obj$par[(n0cum[i] + 1):(n0cum[(i + 1)])]
    gj <- P.posterior[[i]] %*% bj
    F.posterior[, i] <- gj
  }
  B.posterior <- F.posterior %*% solve(A0.mode)
  AR.coefs.posterior <- t(B.posterior[1:(m * p), ])
  dim(AR.coefs.posterior) <- c(m, m, p)
  AR.coefs.posterior <- aperm(AR.coefs.posterior, c(2, 1, 3))
  structural.innovations <- Y1 %*% A0.mode - X1 %*% F.posterior
  if (nexog == 0) {
    exog.coefs <- NA
  }
  else {
    exog.coefs <- B.posterior[((m * p) + 2):nrow(B.posterior), 
                              ]
  }
  gc()
  gc()
  output <- list(XX = XX, XY = XY, YY = YY, y = Y, Y = Y1, 
                 X = X1, structural.innovations = structural.innovations, 
                 Ui = Ui, Hpinv.tilde = Hpinv.tilde, H0inv.tilde = H0inv.tilde, 
                 Pi.tilde = Pi.tilde, Hpinv.posterior = Hpinv.posterior, 
                 P.posterior = P.posterior, H0inv.posterior = H0inv.posterior, 
                 A0.mode = A0.mode, F.posterior = F.posterior, B.posterior = B.posterior, 
                 ar.coefs = AR.coefs.posterior, intercept = B.posterior[(m * 
                                                                           p + 1), ], exog.coefs = exog.coefs, prior = c(lambda0, 
                                                                                                                         lambda1, lambda3, lambda4, lambda5, mu5, mu6), df = capT, 
                 n0 = n0, ident = ident, b = max.obj$par)
  class(output) <- c("BSVAR")
  attr(output, "eqnames") <- colnames(Y)
  return(output)
}


#szbvar
szbvar <- function (Y, p, z = NULL, lambda0, lambda1, lambda3, lambda4, 
          lambda5, mu5, mu6, nu = ncol(Y) + 1, qm = 4, prior = 0, posterior.fit = FALSE) 
{
  sanity.check.bvar(list(Y = Y, p = p, z = z, lambda0 = lambda0, 
                         lambda1 = lambda1, lambda3 = lambda3, lambda4 = lambda4, 
                         lambda5 = lambda5, mu5 = mu5, mu6 = mu6, qm = qm, prior = prior, 
                         posterior.fit = posterior.fit))
  n <- nrow(Y)
  m <- ncol(Y)
  if (is.null(z)) {
    num.exog <- 0
  }
  else {
    num.exog <- ncol(z)
    z <- as.matrix(z)
  }
  ncoef <- (m * p) + 1
  ndum <- m + 1
  capT <- n - p + ndum
  Ts <- n - p
  dat <- as.matrix(Y)
  if (p == 1) {
    datint <- as.vector(dat[1, ])
  }
  else {
    datint <- as.vector(apply(as.matrix(dat[1:p, ]), 2, mean))
  }
  X <- matrix(0, nrow = capT, ncol = ncoef)
  Y <- matrix(0, nrow = capT, ncol = m)
  const <- matrix(1, nrow = capT)
  const[1:m, ] <- 0
  X[, ncoef] <- const
  for (i in 1:m) {
    Y[ndum, i] <- datint[i]
    Y[i, i] <- datint[i]
    for (j in 1:p) {
      X[ndum, m * (j - 1) + i] <- datint[i]
      X[i, m * (j - 1) + i] <- datint[i]
    }
  }
  for (i in 1:p) {
    X[(ndum + 1):capT, (m * (i - 1) + 1):(m * i)] <- matrix(dat[(p + 
                                                                   1 - i):(n - i), ], ncol = m)
  }
  if (is.null(z) == F) {
    pad.z <- matrix(0, nrow = capT, ncol = ncol(z))
    pad.z[(ndum + 1):capT, ] <- matrix(z[(p + 1):n, ], ncol = ncol(z))
    X <- cbind(X, pad.z)
  }
  Y[(ndum + 1):capT, ] <- matrix(dat[(p + 1):n, ], ncol = m)
  X[1:m, ] <- mu5 * X[1:m, ]
  Y[1:m, ] <- mu5 * Y[1:m, ]
  X[ndum, ] <- mu6 * X[ndum, ]
  Y[ndum, ] <- mu6 * Y[ndum, ]
  ld <- seq(1:p)^-lambda3
  if (qm == 12) {
    j <- ceiling(p/3)^-lambda3
    b <- 0
    if (p > 1) {
      b <- (log(1) - log(j))/(1 - p)
    }
    a <- exp(-b)
    ld <- a * exp(b * seq(1:p))
  }
  s2 <- matrix(0, nrow = m, ncol = 1)
  for (i in 1:m) {
    s2[i, 1] <- ar.ols(Y[, i], aic = FALSE, order.max = p, 
                       intercept = TRUE, demean = FALSE)$var.pred
  }
  S0 <- diag(m)
  diag(S0) <- s2/(lambda0^2)
  prior.intercept <- 0
  if (lambda4 > 0) {
    prior.intercept <- 1/(lambda0 * lambda4)^2
  }
  prior.exog <- 0
  if (lambda5 > 0) {
    prior.exog <- 1/(lambda0 * lambda5)^2
  }
  if (num.exog == 0) {
    H0 <- diag(c(kronecker((1/(ld * lambda0 * lambda1))^2, 
                           s2), prior.intercept), nrow = ncoef, ncol = ncoef)
  }
  else {
    H0 <- diag(c(kronecker((1/(ld * lambda0 * lambda1))^2, 
                           s2), prior.intercept, rep(prior.exog, num.exog)), 
               nrow = (ncoef + num.exog), ncol = (ncoef + num.exog))
  }
  if (prior == 1) {
    S0 <- 0 * S0
  }
  if (prior == 2) {
    H0 <- 0 * H0
    S0 <- 0 * S0
    nu <- 0
  }
  XX <- crossprod(X)
  hstar1 <- H0 + XX
  Bh <- solve((hstar1), (crossprod(X, Y) + H0[, 1:m]))
  St <- (S0 + crossprod(Y) + H0[1:m, 1:m] - t(Bh) %*% (hstar1) %*% 
           (Bh))
  Sh <- St/(Ts + nu - m - 1)
  hstarinv <- solve(hstar1)
  vcv.Bh <- kronecker(Sh, hstarinv)
  u <- (Y - X %*% (Bh))
  Sh1 <- crossprod(u)/capT
  intercept <- Bh[ncoef, ]
  ar.coefs <- t(Bh[1:(ncoef - 1), ])
  dim(ar.coefs) <- c(m, m, p)
  ar.coefs <- aperm(ar.coefs, c(2, 1, 3))
  if (is.null(z)) {
    exog.coefs <- NA
  }
  else {
    exog.coefs <- Bh[(ncoef + 1):nrow(Bh), ]
  }
  marg.llf <- NA
  marg.post <- NA
  coef.post <- NA
  pfit <- list(capT = capT, m = m, ncoef = ncoef, num.exog = num.exog, 
               nu = nu, H0 = H0, S0 = S0, Y = Y, X = X, hstar1 = hstar1, 
               Sh = Sh, u = u, Bh = Bh, Sh1 = Sh1)
  output <- list(intercept = intercept, ar.coefs = ar.coefs, 
                 exog.coefs = exog.coefs, Bhat = Bh, vcv = Sh1, vcv.Bh = vcv.Bh, 
                 mean.S = Sh, St = St, hstar = (H0 + XX), hstarinv = hstarinv, 
                 H0 = H0, S0 = S0, residuals = u, X = X, Y = Y, y = dat, 
                 z = z, p = p, num.exog = num.exog, qm = qm, prior.type = prior, 
                 prior = c(lambda0, lambda1, lambda3, lambda4, lambda5, 
                           mu5, mu6, nu), pfit = pfit, marg.llf = marg.llf, 
                 marg.post = marg.post, coef.post = coef.post)
  class(output) <- c("BVAR")
  attr(output, "eqnames") <- colnames(dat)
  if (posterior.fit == T) {
    tmp <- posterior.fit.BVAR(output)
    output$marg.llf <- tmp$data.marg.llf
    output$marg.post <- tmp$data.marg.post
    output$coef.post <- tmp$coef.post
  }
  return(output)
}


#uc.forecast
uc.forecast <- function (varobj, nsteps, burnin, gibbs, exog = NULL) 
{
  if (inherits(varobj, "VAR")) {
    stop("Not implemented for VAR models!\nUse a BVAR with a flat-flat prior if you want this case.\n")
  }
  if (inherits(varobj, "BVAR")) {
    output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs, 
                              exog)
    attr(output, "class") <- c("forecast.VAR")
    return(output)
  }
  if (inherits(varobj, "BSVAR")) {
    stop("Not yet implemented for BSVAR models!\n")
  }
}


#var.lag.specification
var.lag.specification <- function (y, lagmax = 20) 
{
  results <- matrix(0, nrow = lagmax, ncol = 4)
  ldets <- matrix(0, nrow = lagmax, ncol = 4)
  for (p in 1:lagmax) {
    x <- y[(lagmax - p + 1):nrow(y), ]
    tmp <- ar.ols(x, order.max = p, aic = F, intercept = T)
    no.param <- length(tmp$ar) + length(tmp$x.intercept)
    obs <- nrow(y) - lagmax
    ldets[lagmax - p + 1, 2] <- log(det(tmp$var.pred))
    ldets[lagmax - p + 1, 1] <- p
    aic <- ldets[lagmax - p + 1, 2] + no.param * 2/obs
    bic <- ldets[lagmax - p + 1, 2] + no.param * log(obs)/obs
    hq <- ldets[lagmax - p + 1, 2] + 2 * log(log(obs)) * 
      (no.param/obs)
    results[p, 1] <- p
    results[p, 2] <- aic
    results[p, 3] <- bic
    results[p, 4] <- hq
  }
  colnames(results) <- c("Lags", "AIC", "BIC", "HQ")
  for (i in 1:(lagmax - 1)) {
    mcorr <- ncol(y) * ldets[i, 1] + 1
    chi <- (obs - mcorr) * (ldets[i + 1, 2] - ldets[i, 2])
    df <- ncol(y)^2
    chi.p <- 1 - pchisq(chi, df)
    ldets[i, 3] <- chi
    ldets[i, 4] <- chi.p
  }
  colnames(ldets) <- c("Lags", "Log-Det", "Chi^2", "p-value")
  output <- list(ldets = ldets, results = results)
  class(output) <- c("var.lag.specification")
  return(output)
}
