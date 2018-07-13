#' Principal component analysis for Zero-inflated Bivariate Poisson variables
#'
#' This function preforms principal component analysis for pair-wise 
#' Zero-inflated Bivariate Poisson distributed variables.
#' @param data a data frame, providing the data for the principal components analysis
#' @param scale. a logical value, indicating whether the variables should be scaled to have unit variance before the analysis takes place. The default is FALSE.
#' @param maxit a scalar, number of maximum iterations
#' @param pres a scalar, the precision of the EM algorithm
#' @export

zibppca <- function(data, scale. = TRUE, maxit = 30, pres = 1e-06){

if(!is.data.frame(data)){
  print('Error: Data are not a data.frame')
  break
}

p <- dim(data)[2]
n <- dim(data)[1] 

print(paste0('Analyzing data with ',p,' variables and ',n,' observations'))

#check if all numbers are integers
cov.mat <- matrix(0,p,p) 
colnames(cov.mat) <- colnames(data)
rownames(cov.mat) <-  colnames(data)

#Run through all p combinations of variables 
count <- 0 
for(i in 1:p){
  for(j in 1:i){
    #Fit the bivariate Poisson with discrete excess zero distribution on diagonal at 0
    fit <- lm.dibp(x~1,y~1,
                   data = data.frame(y = data[,i],x = data[,j],z1 = rep(1,n),z2 = rep(1,n)),
                   distribution = "discrete", jmax = 0, 
                   maxit = maxit, verbose = FALSE, pres = pres)
    #Extract the lambda3 parameter on log scale, and transform
    cov.mat[i,j] <- exp(coefficients(fit)[3])
  }
}  

cov.mat.tmp <- cov.mat
diag(cov.mat.tmp) <- 0
cov.mat <- cov.mat + t(cov.mat.tmp)

if(!isSymmetric(cov.mat)){
  print('Error: Covariance matrix is not symmetric')
  break
}

#Convert to correlation matrix
if(scale.==TRUE){
  cor.mat <- cov2cor(cov.mat)
} else {
  cor.mat <- cov.mat
}

#Perform eigendecompositon (equivalent to svd) of the correlation matrix 
poisson.svd <- svd(cor.mat)
rownames(poisson.svd$u) <- colnames(data)

zibppca <- list()
zibppca$loadings <- poisson.svd$u
zibppca$eigenvalues <- poisson.svd$d^2
zibppca$scores <- data.matrix(data) %*% poisson.svd$u
return(zibppca)
}

lm.dibp <- function (l1, l2, l1l2 = NULL, l3 = ~1, data, common.intercept = FALSE, 
          zeroL3 = FALSE, distribution = "discrete", jmax = 2, maxit = 300, 
          pres = 1e-08, verbose = getOption("verbose")) 
{
  options(warn = -1)
  templist <- list(l1 = l1, l2 = l2, l1l2 = l1l2, l3 = l3, 
                   data = substitute(data), common.intercept = common.intercept, 
                   zeroL3 = zeroL3, distribution = distribution, jmax = jmax, 
                   maxit = maxit, pres = pres, verbose = verbose)
  tempcall <- as.call(c(expression(lm.dibp), templist))
  rm(templist)
  if (common.intercept) {
    formula1.terms <- "1"
  }
  else {
    formula1.terms <- "internal.data1$noncommon"
  }
  namex <- as.character(l1[2])
  namey <- as.character(l2[2])
  x <- data[, names(data) == namex]
  y <- data[, names(data) == namey]
  n <- length(x)
  lengthprintvec <- 1
  maxy <- max(c(x, y))
  dist <- distribution
  if (charmatch(dist, "poisson", nomatch = 0) == 1) {
    distribution <- 2
  }
  else if (charmatch(dist, "geometric", nomatch = 0) == 1) {
    distribution <- 3
  }
  else if (charmatch(dist, "discrete", nomatch = 0) == 1) {
    distribution <- 1
  }
  if (distribution == 1) {
    dilabel <- paste("Inflation Distribution: Discrete with J=", 
                     jmax)
    if (jmax == 0) {
      theta <- 0
    }
    else {
      theta <- 1:jmax * 0 + 1/(jmax + 1)
    }
    di.f <- function(x, theta) {
      JMAX <- length(theta)
      if (x > JMAX) {
        res <- 0
      }
      else if (x == 0) {
        res <- 1 - sum(theta)
      }
      else {
        res <- theta[x]
      }
      res
    }
  }
  else if (distribution == 2) {
    dilabel <- "Inflation Distribution: Poisson"
    theta <- 1
    di.f <- function(x, theta) {
      if (theta > 0) {
        res <- dpois(x, theta)
      }
      else {
        if (x == 0) {
          res <- 1
        }
        else {
          res <- 1e-12
        }
      }
    }
  }
  else if (distribution == 3) {
    dilabel <- "Inflation Distribution: Geometric"
    theta <- 0.5
    di.f <- function(x, theta) {
      if (theta > 0) {
        if (theta == 1) {
          theta <- 0.9999999
        }
        res <- dgeom(x, theta)
      }
      else if (theta == 1) {
        if (x == 0) {
          res <- 1
        }
        else {
          res <- 1e-12
        }
      }
      else {
        res <- 1e-12
      }
    }
  }
  else {
    stop(paste(distribution, "Not known distribution.", 
               sep = ": "))
  }
  internal.data1 <- data.frame(y1y2 = c(x, y))
  internal.data2 <- data.frame(y3 = rep(0, n))
  p <- length(as.data.frame(data))
  data1 <- rbind(data, data)
  names(data1) <- names(data)
  data1 <- data1[, names(data1) != namex]
  data1 <- data1[, names(data1) != namey]
  if (as.character(l1[3]) == ".") {
    l1 <- formula(paste(as.character(l1[2]), paste(names(data1), 
                                                   "", collapse = "+", sep = ""), sep = "~"))
  }
  if (as.character(l2[3]) == ".") {
    l2 <- formula(paste(as.character(l2[2]), paste(names(data1), 
                                                   "", collapse = "+", sep = ""), sep = "~"))
  }
  if (as.character(l3[2]) == ".") {
    l3 <- formula(paste("", paste(names(data1), "", collapse = "+", 
                                  sep = ""), sep = "~"))
  }
  formula2 <- formula(paste("internal.data2$y3~", as.character(l3[2]), 
                            sep = ""))
  internal.data1$noncommon <- as.factor(c(1:n * 0, 1:n * 0 + 
                                            1))
  contrasts(internal.data1$noncommon) <- contr.treatment(2, 
                                                         base = 1)
  internal.data1$indct1 <- c(1:n * 0 + 1, 1:n * 0)
  internal.data1$indct2 <- c(1:n * 0, 1:n * 0 + 1)
  if (!zeroL3) {
    data2 <- data1[1:n, ]
    names(data2) <- names(data1)
  }
  if (!is.null(l1l2)) {
    formula1.terms <- paste(formula1.terms, as.character(l1l2[2]), 
                            sep = "+")
  }
  templ1 <- labels(terms(l1))
  if (length(templ1) > 0) {
    for (k1 in 1:length(templ1)) {
      if (!is.null(l1l2)) {
        checkvar1 <- sum(labels(terms(l1l2)) == templ1[k1]) == 
          1
      }
      else {
        checkvar1 <- FALSE
      }
      checkvar2 <- sum(labels(terms(l2)) == templ1[k1]) == 
        1
      if (checkvar1 & checkvar2) {
        formula1.terms <- paste(formula1.terms, paste("internal.data1$noncommon*", 
                                                      templ1[k1], sep = ""), sep = "+")
      }
      else {
        formula1.terms <- paste(formula1.terms, paste("+I(internal.data1$indct1*", 
                                                      templ1[k1], sep = ""), sep = "")
        formula1.terms <- paste(formula1.terms, ")", 
                                sep = "")
      }
    }
  }
  templ2 <- labels(terms(l2))
  if (length(templ2) > 0) {
    for (k1 in 1:length(templ2)) {
      if (!is.null(l1l2)) {
        checkvar1 <- (sum(labels(terms(l1l2)) == templ2[k1]) + 
                        sum(labels(terms(l1)) == templ2[k1])) != 2
      }
      else {
        checkvar1 <- TRUE
      }
      if (checkvar1) {
        formula1.terms <- paste(formula1.terms, paste("+I(internal.data1$indct2*", 
                                                      templ2[k1], sep = ""), sep = "")
        formula1.terms <- paste(formula1.terms, ")", 
                                sep = "")
      }
    }
  }
  rm(templ1)
  rm(templ2)
  rm(Checkvar1)
  rm(Checkvar2)
  formula1 <- formula(paste("internal.data1$y1y2~", formula1.terms, 
                            sep = ""))
  tmpform1 <- as.character(formula1[3])
  newformula <- formula1
  while (regexpr("c\\(", tmpform1) != -1) {
    temppos1 <- regexpr("c\\(", tmpform1)[1]
    tempfor <- substring(tmpform1, first = temppos1 + 2)
    temppos2 <- regexpr("\\)", tempfor)[1]
    tempvar <- substring(tempfor, first = 1, last = temppos2 - 
                           1)
    temppos3 <- regexpr(", ", tempvar)[1]
    tempname1 <- substring(tempfor, first = 1, last = temppos3 - 
                             1)
    tempname2 <- substring(tempfor, first = temppos3 + 2, 
                           last = temppos2 - 1)
    tempname2 <- sub("\\)", "", tempname2)
    tempvar1 <- data[, names(data) == tempname1]
    tempvar2 <- data[, names(data) == tempname2]
    data1$newvar1 <- c(tempvar1, tempvar2)
    if (is.factor(tempvar1) & is.factor(tempvar2)) {
      data1$newvar1 <- as.factor(data1$newvar1)
      if (all(levels(tempvar1) == levels(tempvar2))) {
        attributes(data1$newvar1) <- attributes(tempvar1)
      }
    }
    tempvar <- sub(", ", "..", tempvar)
    names(data1)[names(data1) == "newvar1"] <- tempvar
    newformula <- sub("c\\(", "", tmpform1)
    newformula <- sub("\\)", "", newformula)
    newformula <- sub(", ", "..", newformula)
    tmpform1 <- newformula
    formula1 <- formula(paste("internal.data1$y1y2~", newformula, 
                              sep = ""))
  }
  rm(temppos1)
  rm(temppos2)
  rm(temppos3)
  rm(tmpform1)
  rm(tempfor)
  rm(tempvar)
  rm(tempvar1)
  rm(tempvar2)
  rm(tempname1)
  rm(tempname2)
  prob <- 0.2
  s <- rep(0, n)
  vi <- 1:n * 0
  v1 <- 1 - c(vi, vi)
  like <- 1:n * 0
  zero <- (x == 0) | (y == 0)
  if (zeroL3) {
    lambda3 <- rep(0, n)
  }
  else {
    lambda3 <- rep(max(0.1, cov(x, y, use = "complete.obs")), 
                   n)
  }
  internal.data1$v1 <- 1 - c(vi, vi)
  lambda <- glm(formula1, family = poisson, data = data1, 
                weights = internal.data1$v1, maxit = 100)$fitted
  lambda1 <- lambda[1:n]
  lambda2 <- lambda[(n + 1):(2 * n)]
  difllike <- 100
  loglike0 <- 1000
  i <- 0
  ii <- 0
  if (zeroL3) {
    loglike <- rep(0, maxit)
    lambda3 <- 1:n * 0
    while ((difllike > pres) && (i <= maxit)) {
      i <- i + 1
      for (j in 1:n) {
        if (zero[j]) {
          s[j] <- 0
          if (x[j] == y[j]) {
            density.di <- di.f(0, theta)
            like[j] <- log((1 - prob) * exp(-lambda1[j] - 
                                              lambda2[j]) + prob * density.di)
            vi[j] <- prob * density.di * exp(-like[j])
          }
          else {
            like[j] <- log(1 - prob) + log(dpois(x[j], 
                                                 lambda1[j])) + log(dpois(y[j], lambda2[j]))
            vi[j] <- 0
          }
        }
        else {
          if (x[j] == y[j]) {
            density.di <- di.f(x[j], theta)
            like[j] <- log((1 - prob) * dpois(x[j], 
                                              lambda1[j]) * dpois(y[j], lambda2[j]) + 
                             prob * density.di)
            vi[j] <- prob * density.di * exp(-like[j])
          }
          else {
            vi[j] <- 0
            like[j] <- log(1 - prob) + log(dpois(x[j], 
                                                 lambda1[j]) * dpois(y[j], lambda2[j]))
          }
        }
      }
      x1 <- x
      x2 <- y
      loglike[i] <- sum(like)
      difllike <- abs((loglike0 - loglike[i])/loglike0)
      loglike0 <- loglike[i]
      prob <- sum(vi)/n
      if (distribution == 1) {
        if (jmax == 0) {
          theta <- 0
        }
        else {
          for (ii in 1:jmax) {
            temp <- as.numeric(((x == ii) & (y == ii)))
            theta[ii] <- sum(temp * vi)/sum(vi)
          }
        }
      }
      else if (distribution == 2) {
        theta <- sum(vi * x)/sum(vi)
      }
      else if (distribution == 3) {
        theta <- sum(vi)/(sum(vi * x) + sum(vi))
      }
      internal.data1$v1 <- 1 - c(vi, vi)
      internal.data1$v1[(internal.data1$v1 < 0) & (internal.data1$v1 > 
                                                     -1e-10)] <- 0
      x1[(x1 < 0) & (x1 > -1e-10)] <- 0
      x2[(x2 < 0) & (x2 > -1e-10)] <- 0
      internal.data1$y1y2 <- c(x1, x2)
      m <- glm(formula1, family = poisson, data = data1, 
               weights = internal.data1$v1, maxit = 100)
      p3 <- length(m$coef)
      beta <- m$coef
      names(beta) <- newnamesbeta(beta)
      betaparameters <- splitbeta(beta)
      lambda <- fitted(m)
      lambda1 <- lambda[1:n]
      lambda2 <- lambda[(n + 1):(2 * n)]
      if (verbose) {
        printvec <- c(i, beta, 100 * prob, theta, loglike[i], 
                      difllike)
        names(printvec) <- c("iter", names(beta), "Mix.p(%)", 
                             paste("theta", 1:length(theta), sep = ""), 
                             "loglike", "Rel.Dif.loglike")
      }
      else {
        printvec <- c(i, 100 * prob, theta, loglike[i], 
                      difllike)
        names(printvec) <- c("iter", "Mix.p(%)", paste("theta", 
                                                       1:length(theta), sep = ""), "loglike", "Rel.Dif.loglike")
      }
      lengthprintvec <- length(printvec)
    }
    if ((distribution == 1) && (jmax == 0)) {
      noparams <- m$rank + 1
    }
    else {
      noparams <- m$rank + length(theta) + 1
    }
    AIC <- -2 * loglike[i] + noparams * 2
    BIC <- -2 * loglike[i] + noparams * log(2 * n)
    x.mean <- x
    x.mean[x == 0] <- 1e-12
    y.mean <- y
    y.mean[y == 0] <- 1e-12
    AIC.sat <- sum(log(dpois(x, x.mean)) + log(dpois(y, 
                                                     y.mean)))
    BIC.sat <- -2 * AIC.sat + (2 * n) * log(2 * n)
    AIC.sat <- -2 * AIC.sat + (2 * n) * 2
    AICtotal <- c(AIC.sat, AIC)
    BICtotal <- c(BIC.sat, BIC)
    names(AICtotal) <- c("Saturated", "DblPois")
    names(BICtotal) <- c("Saturated", "DblPois")
    allbeta <- c(betaparameters$beta1, betaparameters$beta2)
    names(allbeta) <- c(paste("(l1):", names(betaparameters$beta1), 
                              sep = ""), paste("(l2):", names(betaparameters$beta2), 
                                               sep = ""))
    allparameters <- c(allbeta, prob, theta)
    if (distribution == 1) {
      names(allparameters) <- c(names(allbeta), "p", paste("theta", 
                                                           1:length(theta), sep = ""))
    }
    else {
      names(allparameters) <- c(names(allbeta), "p", "theta")
    }
    fittedval1 <- (1 - prob) * m$fitted[1:n]
    fittedval2 <- (1 - prob) * m$fitted[(n + 1):(2 * n)]
    meandiag <- 0
    if ((distribution == 1) && (jmax > 0)) {
      meandiag <- sum(theta[1:jmax] * 1:jmax)
    }
    else if (distribution == 2) {
      meandiag <- theta
    }
    else if (distribution == 3) {
      meandiag <- (1 - theta)/theta
    }
    fittedval1[x == y] <- prob * meandiag + fittedval1[x == 
                                                         y]
    fittedval2[x == y] <- prob * meandiag + fittedval2[x == 
                                                         y]
    result <- list(coefficients = allparameters, fitted.values = data.frame(x = fittedval1, 
                                                                            y = fittedval2), residuals = data.frame(x = x - 
                                                                                                                      fittedval1, y = y - fittedval2), beta1 = betaparameters$beta1, 
                   beta2 = betaparameters$beta2, p = prob, theta = theta, 
                   diagonal.distribution = dilabel, lambda1 = m$fitted[1:n], 
                   lambda2 = m$fitted[(n + 1):(2 * n)], loglikelihood = loglike[1:i], 
                   parameters = noparams, AIC = AICtotal, BIC = BICtotal, 
                   iterations = i, call = tempcall)
  }
  else {
    loglike <- rep(0, maxit)
    while ((difllike > pres) && (i <= maxit)) {
      i <- i + 1
      for (j in 1:n) {
        if (zero[j]) {
          s[j] <- 0
          if (x[j] == y[j]) {
            density.di <- di.f(0, theta)
            like[j] <- log((1 - prob) * exp(-lambda1[j] - 
                                              lambda2[j] - lambda3[j]) + prob * density.di)
            vi[j] <- prob * density.di * exp(-like[j])
          }
          else {
            like[j] <- log(1 - prob) - lambda3[j] + 
              log(dpois(x[j], lambda1[j])) + log(dpois(y[j], 
                                                       lambda2[j]))
            vi[j] <- 0
          }
        }
        else {
          lbp1 <- pbivpois(x[j] - 1, y[j] - 1, lambda = c(lambda1[j], 
                                                          lambda2[j], lambda3[j]), log = TRUE)
          lbp2 <- pbivpois(x[j], y[j], lambda = c(lambda1[j], 
                                                  lambda2[j], lambda3[j]), log = TRUE)
          s[j] <- exp(log(lambda3[j]) + lbp1 - lbp2)
          if (x[j] == y[j]) {
            density.di <- di.f(x[j], theta)
            like[j] <- log((1 - prob) * exp(lbp2) + 
                             prob * density.di)
            vi[j] <- prob * density.di * exp(-like[j])
          }
          else {
            vi[j] <- 0
            like[j] <- log(1 - prob) + lbp2
          }
        }
      }
      x1 <- x - s
      x2 <- y - s
      loglike[i] <- sum(like)
      difllike <- abs((loglike0 - loglike[i])/loglike0)
      loglike0 <- loglike[i]
      prob <- sum(vi)/n
      if (distribution == 1) {
        if (jmax == 0) {
          theta <- 0
        }
        else {
          for (ii in 1:jmax) {
            temp <- as.numeric(((x == ii) & (y == ii)))
            theta[ii] <- sum(temp * vi)/sum(vi)
          }
        }
      }
      else if (distribution == 2) {
        theta <- sum(vi * x)/sum(vi)
      }
      else if (distribution == 3) {
        theta <- sum(vi)/(sum(vi * x) + sum(vi))
      }
      internal.data2$v1 <- 1 - vi
      internal.data2$v1[(internal.data2$v1 < 0) & (internal.data2$v1 > 
                                                     -1e-10)] <- 0
      internal.data2$y3 <- s
      m0 <- glm(formula2, family = poisson, data = data2, 
                weights = internal.data2$v1, maxit = 100)
      beta3 <- m0$coef
      lambda3 <- m0$fitted
      internal.data1$v1 <- 1 - c(vi, vi)
      internal.data1$v1[(internal.data1$v1 < 0) & (internal.data1$v1 > 
                                                     -1e-10)] <- 0
      x1[(x1 < 0) & (x1 > -1e-10)] <- 0
      x2[(x2 < 0) & (x2 > -1e-10)] <- 0
      internal.data1$y1y2 <- c(x1, x2)
      m <- glm(formula1, family = poisson, data = data1, 
               weights = internal.data1$v1, maxit = 100)
      p3 <- length(m$coef)
      beta <- m$coef
      names(beta) <- newnamesbeta(beta)
      lambda <- fitted(m)
      lambda1 <- lambda[1:n]
      lambda2 <- lambda[(n + 1):(2 * n)]
      if (verbose) {
        printvec <- c(i, beta, beta3, 100 * prob, theta, 
                      loglike[i], difllike)
        names(printvec) <- c("iter", names(beta), paste("l3_", 
                                                        names(beta3), sep = ""), "Mix.p(%)", paste("theta", 
                                                                                                   1:length(theta), sep = ""), "loglike", "Rel.Dif.loglike")
      }
      else {
        printvec <- c(i, 100 * prob, theta, loglike[i], 
                      difllike)
        names(printvec) <- c("iter", "Mix.p(%)", paste("theta", 
                                                       1:length(theta), sep = ""), "loglike", "Rel.Dif.loglike")
      }
      lengthprintvec <- length(printvec)
    }
    if ((distribution == 1) && (jmax == 0)) {
      noparams <- m$rank + m0$rank + 1
    }
    else {
      noparams <- m$rank + m0$rank + length(theta) + 1
    }
    AIC <- -2 * loglike[i] + noparams * 2
    BIC <- -2 * loglike[i] + noparams * log(2 * n)
    x.mean <- x
    x.mean[x == 0] <- 1e-12
    y.mean <- y
    y.mean[y == 0] <- 1e-12
    AIC.sat <- sum(log(dpois(x, x.mean)) + log(dpois(y, 
                                                     y.mean)))
    BIC.sat <- -2 * AIC.sat + (2 * n) * log(2 * n)
    AIC.sat <- -2 * AIC.sat + (2 * n) * 2
    AICtotal <- c(AIC.sat, AIC)
    BICtotal <- c(BIC.sat, BIC)
    names(AICtotal) <- c("Saturated", "BivPois")
    names(BICtotal) <- c("Saturated", "BivPois")
    betaparameters <- splitbeta(beta)
    allbeta <- c(betaparameters$beta1, betaparameters$beta2, 
                 beta3)
    names(allbeta) <- c(paste("(l1):", names(betaparameters$beta1), 
                              sep = ""), paste("(l2):", names(betaparameters$beta2), 
                                               sep = ""), paste("(l3):", names(beta3), sep = ""))
    allparameters <- c(allbeta, prob, theta)
    if (distribution == 1) {
      names(allparameters) <- c(names(allbeta), "p", paste("theta", 
                                                           1:length(theta), sep = ""))
    }
    else {
      names(allparameters) <- c(names(allbeta), "p", "theta")
    }
    fittedval1 <- (1 - prob) * (m$fitted[1:n] + lambda3)
    fittedval2 <- (1 - prob) * (m$fitted[(n + 1):(2 * n)] + 
                                  lambda3)
    meandiag <- 0
    if ((distribution == 1) && (jmax > 0)) {
      meandiag <- sum(theta[1:jmax] * 1:jmax)
    }
    else if (distribution == 2) {
      meandiag <- theta
    }
    else if (distribution == 3) {
      meandiag <- (1 - theta)/theta
    }
    fittedval1[x == y] <- prob * meandiag + fittedval1[x == 
                                                         y]
    fittedval2[x == y] <- prob * meandiag + fittedval2[x == 
                                                         y]
    result <- list(coefficients = allparameters, fitted.values = data.frame(x = fittedval1, 
                                                                            y = fittedval2), residuals = data.frame(x = x - 
                                                                                                                      fittedval1, y = y - fittedval2), beta1 = betaparameters$beta1, 
                   beta2 = betaparameters$beta2, beta3 = beta3, p = prob, 
                   theta = theta, diagonal.distribution = dilabel, 
                   lambda1 = m$fitted[1:n], lambda2 = m$fitted[(n + 
                                                                  1):(2 * n)], lambda3 = lambda3, loglikelihood = loglike[1:i], 
                   parameters = noparams, AIC = AICtotal, BIC = BICtotal, 
                   iterations = i, call = tempcall)
  }
  options(warn = 0)
  class(result) <- c("lm.dibp", "lm")
  result
}

pbivpois <- function (x, y = NULL, lambda = c(1, 1, 1), log = FALSE) 
{
  if (is.matrix(x)) {
    var1 <- x[, 1]
    var2 <- x[, 2]
  }
  else if (is.vector(x) & is.vector(y)) {
    if (length(x) == length(y)) {
      var1 <- x
      var2 <- y
    }
    else {
      stop("lengths of x and y are not equal")
    }
  }
  else {
    stop("x is not a matrix or x and y are not vectors")
  }
  n <- length(var1)
  logbp <- vector(length = n)
  for (k in 1:n) {
    x0 <- var1[k]
    y0 <- var2[k]
    xymin <- min(x0, y0)
    lambdaratio <- lambda[3]/(lambda[1] * lambda[2])
    i <- 0:xymin
    sums <- -lgamma(var1[k] - i + 1) - lgamma(i + 1) - lgamma(var2[k] - 
                                                                i + 1) + i * log(lambdaratio)
    maxsums <- max(sums)
    sums <- sums - maxsums
    logsummation <- log(sum(exp(sums))) + maxsums
    logbp[k] <- -sum(lambda) + var1[k] * log(lambda[1]) + 
      var2[k] * log(lambda[2]) + logsummation
  }
  if (log) {
    result <- logbp
  }
  else {
    result <- exp(logbp)
  }
  result
}

newnamesbeta <- function (bvec) 
{
  names(bvec) <- sub("\\)", "", names(bvec))
  names(bvec) <- sub("\\(Intercept", "(Intercept)", names(bvec))
  names(bvec)[pmatch("internal.data1$noncommon2", names(bvec))] <- "(l2-l1):(Intercept)"
  names(bvec) <- sub("internal.data1\\$noncommon2:", "(l2-l1):", 
                     names(bvec))
  names(bvec) <- sub("internal.data1\\$noncommon0:", "(l1):", 
                     names(bvec))
  names(bvec) <- sub("internal.data1\\$noncommon1:", "(l2):", 
                     names(bvec))
  names(bvec) <- sub(":internal.data1\\$noncommon2", "(l2-l1):", 
                     names(bvec))
  names(bvec) <- sub(":internal.data1\\$noncommon0", "(l1):", 
                     names(bvec))
  names(bvec) <- sub(":internal.data1\\$noncommon1", "(l2):", 
                     names(bvec))
  names(bvec) <- sub("I\\(internal.data1\\$indct1 \\* ", "(l1):", 
                     names(bvec))
  names(bvec) <- sub("I\\(internal.data1\\$indct2 \\* ", "(l2):", 
                     names(bvec))
  names(bvec)
}

splitbeta <- function (bvec) 
{
  p3 <- length(bvec)
  indx1 <- grep("\\(l1\\):", names(bvec))
  indx2 <- grep("\\(l2\\):", names(bvec))
  indx3 <- grep("\\(l2-l1\\):", names(bvec))
  tempnames <- sub("\\(l2-l1)\\:", "k", names(bvec))
  tempnames <- sub("\\(l2)\\:", "k", tempnames)
  tempnames <- sub("\\(l1)\\:", "k", tempnames)
  indx4 <- tempnames %in% names(bvec)
  beta1 <- c(bvec[indx4], bvec[indx1])
  beta2 <- c(bvec[indx4], bvec[indx3], bvec[indx2])
  indexbeta2 <- c(rep(0, sum(indx4)), rep(1, length(indx3)), 
                  rep(2, length(indx2)))
  names(beta1) <- sub("\\(l1\\):", "", names(beta1))
  names(beta2) <- sub("\\(l2\\):", "", names(beta2))
  names(beta2) <- sub("\\(l2-l1\\):", "", names(beta2))
  beta1 <- beta1[order(names(beta1))]
  indexbeta2 <- indexbeta2[order(names(beta2))]
  beta2 <- beta2[order(names(beta2))]
  ii <- 1:length(beta2)
  ii <- ii[indexbeta2 == 0]
  for (i in ii) {
    beta2[i] <- sum(beta2[names(beta2)[i] == names(beta2)])
  }
  beta2 <- beta2[indexbeta2 %in% c(0, 2)]
  btemp <- list(beta1 = beta1, beta2 = beta2)
  btemp
}