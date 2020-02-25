
# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
HMMsegment <- function(correctOut, param = NULL, autosomes = NULL, maxiter = 50,
    getparam = FALSE, verbose = TRUE) {
  chr <- correctOut$chr

  if (is.null(autosomes)) {
    autosomes <- (chr != "X" & chr != "Y" & chr != "23" & chr != "24" &
      chr != "chrX" & chr != "chrY" & chr != "M" & chr != "MT" & chr != "chrM")
  }

  if (is.null(param)) {
    param <- data.frame(
      strength = 10000000,
      e = 0.9999999,
      mu = quantile(correctOut$copy, na.rm = TRUE,
        prob = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)),
      lambda = 20,
      nu = 2.1,
      kappa = c(0.05, 0.05, 0.7, 0.1, 0.05, 0.05) * 1000,
      m = 0,
      eta = c(5, 5, 50, 5, 5, 5) * 10000,
      gamma = 3,
      S = 0
    )
    param$m <- param$mu
    param$S <-
      ((sd(2 ^ correctOut$copy[autosomes], na.rm = TRUE)/sqrt(nrow(param))) ^ 2)
    rownames(param) <- seq(1, 6)
  }

  if (getparam) {
    return(param)
  }

  output <- manualSegment(correctOut$copy, chr, autosomes, param, maxiter,
    verbose)
  output$segs <- processSegments(output$segs, chr, correctOut$start,
    correctOut$end, correctOut$copy)
  return(output)
}

# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
processSegments <- function(seg, chr, start, end, copy) {
  segment <- data.frame()
  chromosomes <- levels(chr)
  for (i in 1:length(chromosomes)) {
    seg_length = dim(seg[[i]])[1]
    chr_name <- rep(chromosomes[i], seg_length)
    chr_index <- which(chr == chromosomes[i])
    chr_start <- start[chr_index][seg[[i]][, 1]]
    chr_stop <- end[chr_index][seg[[i]][, 2]]
    chr_state <- seg[[i]][, 3]
    chr_median <- rep(0, seg_length)
    for(j in 1:seg_length) {
      chr_median[j] <-
        median(na.rm = TRUE, copy[chr_index][seg[[i]][j, 1]:seg[[i]][j, 2]])
    }
    segment <- rbind(segment, cbind(chr = chr_name,
      start = as.numeric(chr_start), end = chr_stop, state = chr_state,
      median = chr_median))
  }
  segment <- transform(segment, start = as.numeric(as.character(start)),
    end = as.numeric(as.character(end)),
    median = as.numeric(as.character(median)))
  return(segment)
}

# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
manualSegment <- function(copy, chr, autosomes, param, maxiter,
    verbose = TRUE) {

  if (length(copy) != length(chr) || length (copy) != length(autosomes)) {
    stop("Length of inputs do not match for one of: copy, chr, autosomes")
  }

  if (is.null(param$mu) || is.null(param$lambda) || is.null(param$nu) ||
    is.null(param$kappa) || is.null(param$eta) || is.null(param$gamma)) {
    stop(paste("Parameter missing, ensure all parameters exist as columns in",
      "data frame: mu, lambda, nu, kappa, eta, gamma"))
  }

  K <- dim(param)[1]
  N <- length(copy)
  rho <- matrix(0, K, N)
  py <- matrix(0, K, N)            # Local evidence
  mus <- matrix(0, K, maxiter)     # State means
  lambdas <- matrix(0, K, maxiter) # State Variances
  kappa <- param$kappa             # Initial state distribution
  converged <- FALSE               # Flag for convergence
  Z <- rep(0, N)
  Zcounts <- matrix(0, K, K)
  loglik <- rep(0, maxiter)

  # SET UP
  # Set up the chromosome indices and make cell array of chromosome indicies
  chrs <- levels(chr)              # WARNING, gets resorted here...
  chrsI <- vector('list', length(chrs))

  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }

  # INITIALIZATION
  if (verbose) { message("Initialization") }
  i <- 1
  mus[, 1] <- param$mu # First column
  lambdas[, 1] <- param$lambda
  # Recalculate the likelihood
  for (k in 1:K) {
    py[k, ] <- tdistPDF(copy, mus[k, i],lambdas[k, i], param$nu[k]);
  }
  # Initialize transition matrix to the prior
  A <- matrix(0, K, K)
  for (j in 1:K) {
    A[j, ] <- (1 - param$e[1]) / (K - 1)
    A[j, j] <- param$e[1]
  }
  A_prior <- A
  dirPrior <- A * param$strength[1]
  loglik[i] <- -Inf

  while(!converged && (i < maxiter)) {
    if (verbose) { message(paste("EM iteration:", i,
      "Log likelihood:", loglik[i])) }
    i <- i + 1

    # E-step
    ############################################################################
    if (verbose) { message("Expectation") }
    Zcounts <- matrix(0, K, K)
    for (j in 1:length(chrsI)) {
      I <- chrsI[[j]]
      output <- .Call("forward_backward", kappa, A, py[, I],PACKAGE = "CaSpER")
      rho[, I] <- output$rho
      loglik[i] <- loglik[i] + output$loglik
      Zcounts <- Zcounts + t(colSums(aperm(output$xi, c(3, 2, 1))))
    }

    # M-step
    ############################################################################
    # Update the noise hyperparams
    if (verbose) { message("Maximization") }
    mu_i <- mus[, i - 1]
    lambda_i <- lambdas[, i - 1]
    output <- estimateTNoiseParamsMap(copy[autosomes], mu_i, lambda_i, param$nu,
      rho[, autosomes], param$eta, param$m, param$gamma, param$S, param$kappa)
    mus[, i] <- output$mu_N
    lambdas[, i] <- output$lambda_N
    kappa <- output$pi

    # Recalculate the likelihood
    for (k in 1:K) {
      py[k, ] <- tdistPDF(copy, mus[k, i], lambdas[k, i], param$nu[k])
    }

    # Update transition matrix A
    priorA <- 0
    for (k in 1:K) {
      A[k, ] <- Zcounts[k, ] + dirPrior[k, ]
      A[k, ] <- normalize(A[k, ])
      priorA <- priorA + log(dirichletpdf(A_prior[k, ], A[k, ]))
    }

    # Compute log likelihood and check convergence
    priorMu <- c()
    for(k in 1:K) {
      priorMu[k] <- log(dnorm(mus[k, i], param$mu[k], 1))
    }
    loglik[i] <- loglik[i] + priorA + sum(priorMu)
    if (abs(loglik[i] - loglik[i - 1]) < 1e-1 || loglik[i] < loglik[i - 1]) {
      converged = 1
    }
  }

  if (converged) {
    i = i - 1
     # Perform one last round of E-step to get latest responsibilities
     #E-step
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (verbose) { message("Re-calculating latest responsibilties for output") }   
    for (j in 1:length(chrsI)) {
      I <- chrsI[[j]]
      output <- .Call("forward_backward", kappa, A, py[, I],
        PACKAGE = "CaSpER")
      rho[, I] <- output$rho
    }

  }

  if (verbose) {
    message("Optimal parameters found, segmenting and classifying")
  }

  segs <- vector('list', length(chrs))
  for(c in 1:length(chrsI)) {
    I <- chrsI[[c]]
    output <- .Call("viterbi", log(kappa), log(A), log(py[, I]),PACKAGE = "CaSpER")
    Z[I] <- output$path
    segs[[c]] <- output$seg
  }

  mus = mus[, 1:i]
  lambdas = lambdas[, 1:i]
  loglik = loglik[1:i]

  output <- vector('list', 0);
  output$state <- Z
  output$segs <- segs
  output$mus <- mus
  output$lambdas <- lambdas
  output$pi <- kappa
  output$loglik <- loglik
  output$rho <- rho

  return(output)
}


# Normalize a given array to sum to 1
# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
normalize <- function(A) {
  M <- A / (sum(A) + (sum(A) == 0))
  return(M);
}

# Dirichlet probability density function, returns the probability of vector
# x under the Dirichlet distribution with parameter vector alpha
# Author: David Ross http://www.cs.toronto.edu/~dross/code/dirichletpdf.m
dirichletpdf <- function(x, alpha) {
  if (any(x < 0)) {
    return(0);
  }
  if (abs(sum(x) - 1) > 1e-3) {
    stop("Dirichlet PDF: probabilities do not sum to 1")
    return(0);
  }
  p <- exp(lgamma(sum(alpha)) - sum(lgamma(alpha))) * prod(x ^ (alpha - 1))
  return(p)
}

# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
tdistPDF <- function(x, mu, lambda, nu) {
  p <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda / (pi * nu)) ^ (0.5)) *
    (1 + (lambda * (x - mu) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
  p[is.na(p)] <- 1
  return(p)
}

# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
estimateTNoiseParamsMap <- function(
    y, mu, lambda, nu, rho, eta, m, gamma, S, kappa) {

  yr <- t(matrix(y, length(y), length(mu)))        # Vectorize parameter
  u <- (1 + nu) / (((yr - mu) ^ 2) * lambda + nu); # scale parameter

  # Calculate the mean
  mu_N = (rowSums(rho * u * yr, na.rm = TRUE) + (eta * m)) /
    (rowSums(rho * u, na.rm = TRUE) + eta)

  # Calculate the precision
  lambda_N = (rowSums(rho, na.rm = TRUE) + gamma + 1) /
    (rowSums(rho * u * ((yr - mu_N) ^ 2), na.rm = TRUE) +
    (eta * (mu_N - m) ^ 2 + S))

  # Calculate the stationary distribution
  pi <- (rowSums(rho, na.rm = TRUE) + kappa - 1) /
    (length(y) + sum(kappa, na.rm = TRUE) + length(mu))

  out <- vector('list', 0)
  out$mu_N <- mu_N
  out$lambda_N <- lambda_N
  out$pi <- pi
  return(out)
}

