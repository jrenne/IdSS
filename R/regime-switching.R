# Markov-switching Gaussian AR(1) models. No additional dependency.
# Convention: P[i,j] = Pr(s_t=j | s_(t-1)=i), so rows sum to one.
# c is a recursion intercept, not an unconditional regime mean.

ms_check_model <- function(c, phi, sigma, P) {
  K <- length(c)
  if (!.idss_real(c) || !is.null(dim(c)) || !.idss_real(phi) ||
      !is.null(dim(phi)) || !.idss_real(sigma) || !is.null(dim(sigma)) ||
      !.idss_real(P) || K < 1L || length(phi) != K || length(sigma) != K ||
      !is.matrix(P) || !identical(dim(P),c(K,K)) ||
      any(!is.finite(c(c,phi,sigma,P))) || any(sigma <= 0) ||
      any(P < 0) || max(abs(rowSums(P)-1)) > 1e-10)
    stop("Invalid MS-AR parameters or row-stochastic transition matrix")
  invisible(K)
}

#' Stationary probabilities of a finite Markov chain
#'
#' Compute the unique stationary distribution of a row-stochastic finite-state transition matrix.
#'
#' @param P Finite square row-stochastic matrix.
#' @return Numeric column-probability vector q satisfying t(P) %*% q = q and sum(q) = 1; stored as an R vector.
#' @details Irreducibility is sufficient, but not necessary, for uniqueness. Periodicity does not prevent a unique stationary distribution. An error is raised if no unique solution exists. No initial distribution is imposed on a filter by calling this function.
#' @examples
#' P <- matrix(c(.9,.1,.05,.95),2,byrow=TRUE)
#' q <- ms_stationary(P)
#' stopifnot(max(abs(as.vector(t(P) %*% q)-q)) < 1e-12)
#' @family Markov-switching models
#' @export
ms_stationary <- function(P) {
  if (!is.matrix(P) || !.idss_real(P) || !nrow(P) || nrow(P)!=ncol(P) ||
      any(!is.finite(P)) || any(P<0) || max(abs(rowSums(P)-1))>1e-10)
    stop("P must be a finite square row-stochastic matrix")
  K <- nrow(P)
  A <- t(P)-diag(K); A[K,] <- 1
  b <- c(rep(0,K-1L),1)
  q <- tryCatch(as.vector(solve(A,b)),error=function(e)
    stop("P must have a unique stationary distribution"))
  if (min(q) < -1e-10) stop("No admissible unique stationary distribution")
  q <- pmax(q,0); q/sum(q)
}

# The prior is Pr(s_1 | y_0). This is a specified conditional initialization.
#' Hamilton filtering and backward smoothing for a switching AR(1)
#'
#' Evaluate the conditional Gaussian likelihood and state probabilities of a finite-state Markov-switching AR(1), optionally using future observations for backward smoothing.
#'
#' @param c Numeric vector of regime-specific recursion intercepts, not unconditional regime means.
#' @param phi Numeric vector of regime-specific AR(1) coefficients.
#' @param sigma Positive numeric vector of regime innovation standard deviations.
#' @param P Square row-stochastic transition matrix: P[i,j] is the probability of moving from state i to state j. Probability vectors are columns.
#' @param y Nonempty finite numeric observation vector. Missing values are not supported by this routine.
#' @param y0 Finite pre-sample observation on which the likelihood conditions.
#' @param prior Initial probabilities Pr(s[1] | y0). NULL uses uniform probabilities, held fixed for this call.
#' @param smooth Logical scalar. If TRUE, also run the exact finite-state backward smoother.
#' @return A list containing loglik, per-observation increments, n-by-K predicted and filtered probability matrices, y, y0, c, phi, sigma, P and prior. With smooth=TRUE it also contains the n-by-K smoothed matrix and a K-by-K matrix of expected transition_counts over dates 2 through n.
#' @details The observation equation is y[t] = c[s[t]] + phi[s[t]] * y[t-1] + sigma[s[t]] * z[t], where the z innovations are independent standard Gaussians and independent of the autonomous Markov chain. This is a switching-intercept model, not the switching-mean AR model whose density depends on lagged regimes. Computation uses logarithmic normalization. Zero transition probabilities are allowed. Smoothing conditions on the complete observed sample and fixed parameters; it is not available in real time. A stationary chain prior does not make a fixed y0 a draw from the joint stationary law.
#' @examples
#' P <- matrix(c(.9,.1,.2,.8),2,byrow=TRUE)
#' f <- ms_filter(c(1,-1,.2),c(-1,1),c(0,0),c(1,1),P,prior=c(.5,.5))
#' f$filtered
#' f$smoothed
#' stopifnot(abs(sum(f$transition_counts)-2) < 1e-10)
#' @family Markov-switching models
#' @export
#' @importFrom stats dnorm
#' @importFrom utils head
ms_filter <- function(y, c, phi, sigma, P, y0=0, prior=NULL, smooth=TRUE) {
  K <- ms_check_model(c,phi,sigma,P)
  if (!.idss_real(y) || !is.null(dim(y)) || !.idss_flag(smooth))
    stop("y must be a numeric vector and smooth a logical scalar")
  y <- as.numeric(y); n <- length(y)
  if (!n || any(!is.finite(y)) || !.idss_scalar(y0))
    stop("Complete finite observations and a scalar y0 are required")
  if (is.null(prior)) prior <- rep(1/K,K)
  if (!.idss_prob(prior,K)) stop("Invalid initial state probabilities")
  lag <- c(y0,head(y,-1))
  log_density <- vapply(seq_len(K),function(j)
    dnorm(y,c[j]+phi[j]*lag,sigma[j],log=TRUE),numeric(n))
  log_density <- matrix(log_density,n,K)
  predicted <- filtered <- matrix(0,n,K)
  increments <- numeric(n)
  for (t in seq_len(n)) {
    predicted[t,] <- if(t==1L) prior else as.vector(crossprod(P,filtered[t-1,]))
    a <- log(predicted[t,])+log_density[t,]
    shift <- max(a)
    if (!is.finite(shift)) stop("Non-finite predictive log density at observation ",t)
    w <- exp(a-shift); total <- sum(w)
    increments[t] <- shift+log(total)
    filtered[t,] <- w/total
  }
  result <- list(loglik=sum(increments),increments=increments,
                 predicted=predicted,filtered=filtered,y=y,y0=y0,
                 c=c,phi=phi,sigma=sigma,P=P,prior=prior)
  if (smooth) {
    smoothed <- filtered
    transition_counts <- matrix(0,K,K)
    if(n>1L) for(t in (n-1L):1L) {
      ratio <- numeric(K); possible <- predicted[t+1,]>0
      ratio[possible] <- smoothed[t+1,possible]/predicted[t+1,possible]
      pair <- P*outer(filtered[t,],ratio)
      pair <- pair/sum(pair)
      smoothed[t,] <- rowSums(pair)
      transition_counts <- transition_counts+pair
    }
    result$smoothed <- smoothed
    result$transition_counts <- transition_counts
  }
  result
}

#' Simulate a Gaussian Markov-switching AR(1)
#'
#' Simulate observations and their underlying finite-state regime sequence, with optional burn-in.
#'
#' @param c Numeric vector of regime-specific recursion intercepts, not unconditional regime means.
#' @param phi Numeric vector of regime-specific AR(1) coefficients.
#' @param sigma Positive numeric vector of regime innovation standard deviations.
#' @param P Square row-stochastic transition matrix: P[i,j] is the probability of moving from state i to state j. Probability vectors are columns.
#' @param n Positive integer number of retained observations.
#' @param burn Non-negative integer number of observations discarded before retention.
#' @param y0 Finite scalar starting observation before burn-in.
#' @return A list containing the retained numeric y vector, integer state labels from 1 through K, and the observation y0 immediately preceding the retained sample.
#' @details The first simulated regime is drawn from the unique stationary distribution of P. The observation starts from the specified y0. A finite burn-in approximates stationary observation initialization in a stable model; it is not exact. Uniform contraction max(abs(phi)) < 1 is a sufficient finite-variance stability condition, but simulation does not impose it.
#' @examples
#' set.seed(123)
#' P <- matrix(c(.9,.1,.05,.95),2,byrow=TRUE)
#' x <- ms_simulate(100,c(-1,1),c(.3,.3),c(.5,.5),P,burn=200)
#' plot(x$y,type="l",xlab="Time",ylab="Observation")
#' table(x$state)
#' @family Markov-switching models
#' @export
ms_simulate <- function(n, c, phi, sigma, P, burn=1000L, y0=0) {
  K <- ms_check_model(c,phi,sigma,P)
  if(!.idss_integer(n) || !.idss_integer(burn,0L) ||
     n+burn>.Machine$integer.max || !.idss_scalar(y0))
    stop("n and burn must be admissible integers")
  total <- n+burn
  state <- integer(total); y <- numeric(total)
  current <- sample.int(K,1,prob=ms_stationary(P)); previous <- y0
  for(t in seq_len(total)) {
    if(t>1L) current <- sample.int(K,1,prob=P[current,])
    state[t] <- current
    y[t] <- c[current]+phi[current]*previous+sigma[current]*rnorm(1)
    previous <- y[t]
  }
  take <- burn+seq_len(n)
  list(y=y[take],state=state[take],y0=if(burn>0) y[burn] else y0)
}

# Two-state EM/ECM fit: common AR coefficient, common or switching variance.
# Initialization probabilities stay fixed throughout EM, preserving its likelihood.
# Bounds prevent variance collapse, near-unit-root dynamics and empty regimes.
#' Estimate a two-state switching-intercept AR(1) by EM
#'
#' Fit a two-regime Gaussian model with a common AR coefficient and common or regime-specific innovation variance, using constrained multistart EM or ECM.
#'
#' @param y Finite numeric observation vector, with at least 30 entries and positive sample standard deviation.
#' @param ar Logical scalar. TRUE estimates a common AR(1) coefficient, conditions on the first observation and fits the remaining observations. FALSE fixes the AR coefficient at zero and fits every observation.
#' @param switching_variance Logical scalar selecting two innovation variances rather than one common variance.
#' @param zero_mean Logical scalar fixing both intercepts at zero; requires ar=FALSE.
#' @param starts Positive integer number of deterministic starting configurations.
#' @param maxit Positive integer maximum iterations for each start.
#' @param tol Positive relative log-likelihood convergence tolerance.
#' @param sd_floor Positive lower bound on each innovation standard deviation, smaller than the standard deviation of the fitted observations. NULL uses 0.1 times that standard deviation.
#' @param probability_floor Scalar strictly between zero and 0.5; all transition probabilities lie between this bound and its complement.
#' @param phi_bound Positive scalar below one bounding the absolute common AR coefficient.
#' @return The selected converged fit as returned by ms_filter, plus index (positions in the input y), converged, history, iterations, successful_starts, start_loglik, sd_floor, probability_floor, boundary, npar, aic and bic. The boundary flag indicates a fitted active constraint. start_loglik includes all starts, including starts that did not converge; the selected fit is the best converged start.
#' @details Initial state probabilities remain fixed at (0.5,0.5) in every iteration. Smoothed state and pair probabilities provide the EM weights. With switching variances, conditional maximizations update regression and variance blocks. The largest likelihood among converged starts is retained; no global-maximum guarantee or standard errors are provided. The default labels order intercepts increasingly, or variances increasingly for a zero-mean model. The explicit bounds prevent variance collapse and nearly empty regimes; inspect boundary before interpreting an estimate as an interior optimum. AIC/BIC use the same conditional sample and do not make regular chi-squared tests of the number of regimes valid.
#' @examples
#' set.seed(123)
#' P <- matrix(c(.9,.1,.1,.9),2,byrow=TRUE)
#' y <- ms_simulate(180,c(-1,1),c(.25,.25),c(.5,.5),P,burn=300)$y
#' fit <- fit_msar(y,starts=2)
#' fit$c
#' fit$P
#' fit$boundary
#' stopifnot(fit$converged)
#' @family Markov-switching models
#' @export
fit_msar <- function(y, ar=TRUE, switching_variance=FALSE, zero_mean=FALSE,
                     starts=6L, maxit=1000L, tol=1e-8,
                     sd_floor=NULL, probability_floor=.005, phi_bound=.98) {
  if(!.idss_real(y) || !is.null(dim(y)) || !.idss_flag(ar) ||
     !.idss_flag(switching_variance) || !.idss_flag(zero_mean) ||
     !.idss_integer(starts) || !.idss_integer(maxit) ||
     !.idss_scalar(tol) || tol<=0 || !.idss_scalar(probability_floor) ||
     !.idss_scalar(phi_bound)) stop("Invalid fitting inputs or controls")
  y <- as.numeric(y)
  if(length(y)<30 || any(!is.finite(y)) || sd(y)<=0 ||
     starts<1 || starts!=as.integer(starts) || maxit<1 ||
     probability_floor<=0 || probability_floor>=.5 || phi_bound<=0 || phi_bound>=1)
    stop("Invalid data or fitting controls")
  if(zero_mean && ar) stop("The zero-mean teaching specification uses ar=FALSE")
  index <- if(ar) 2:length(y) else seq_along(y)
  obs <- y[index]; n <- length(obs)
  y0 <- if(ar) y[1] else 0
  lag <- if(ar) head(y,-1) else rep(0,n)
  scale <- sd(obs)
  if(is.null(sd_floor)) sd_floor <- .1*scale
  if(!.idss_scalar(sd_floor) || sd_floor<=0 || sd_floor>=scale)
    stop("sd_floor must be positive and below the sample standard deviation")
  prior <- c(.5,.5)
  candidates <- vector("list",starts)
  persistence <- c(.8,.95,.6,.9,.7,.97)
  for(start in seq_len(starts)) {
    p <- persistence[(start-1L)%%length(persistence)+1L]
    p <- min(1-probability_floor,max(probability_floor,p))
    P <- matrix(c(p,1-p,1-p,p),2,2,byrow=TRUE)
    phi <- if(ar) c(0,.3,-.2,.6)[(start-1L)%%4L+1L] else 0
    phi <- min(phi_bound,max(-phi_bound,phi))
    intercept <- if(zero_mean) c(0,0) else
      as.numeric(quantile(obs,c(.2,.8)))*(1-phi)
    sigma <- if(switching_variance) scale*c(.55,1.5) else rep(scale,2)
    sigma <- pmax(sd_floor,sigma)
    current <- ms_filter(obs,intercept,rep(phi,2),sigma,P,y0,prior)
    history <- current$loglik; converged <- FALSE
    for(iter in seq_len(maxit)) {
      weights <- current$smoothed
      occupancy <- colSums(weights)
      if(min(occupancy)<1e-7) break
      counts <- current$transition_counts
      if(min(rowSums(counts))<1e-10) break
      p11 <- counts[1,1]/sum(counts[1,])
      p22 <- counts[2,2]/sum(counts[2,])
      p11 <- min(1-probability_floor,max(probability_floor,p11))
      p22 <- min(1-probability_floor,max(probability_floor,p22))
      new_P <- matrix(c(p11,1-p11,1-p22,p22),2,2,byrow=TRUE)
      new_phi <- phi; new_c <- intercept
      if(!zero_mean) {
        # Weighted within-regime demeaning profiles out the two intercepts.
        w <- sweep(weights,2,sigma^2,"/")
        sw <- colSums(w)
        xb <- colSums(w*lag)/sw; yb <- colSums(w*obs)/sw
        if(ar) {
          numerator <- denominator <- 0
          for(j in 1:2) {
            numerator <- numerator+sum(w[,j]*(lag-xb[j])*(obs-yb[j]))
            denominator <- denominator+sum(w[,j]*(lag-xb[j])^2)
          }
          if(denominator<=0) break
          new_phi <- min(phi_bound,max(-phi_bound,numerator/denominator))
        } else new_phi <- 0
        new_c <- yb-new_phi*xb
      }
      residual <- outer(obs-new_phi*lag,c(1,1))-matrix(new_c,n,2,byrow=TRUE)
      if(switching_variance) new_sigma <- sqrt(colSums(weights*residual^2)/occupancy)
      else new_sigma <- rep(sqrt(sum(weights*residual^2)/n),2)
      new_sigma <- pmax(sd_floor,new_sigma)
      next_fit <- ms_filter(obs,new_c,rep(new_phi,2),new_sigma,new_P,y0,prior)
      gain <- next_fit$loglik-current$loglik
      if(gain < -1e-7*(1+abs(current$loglik))) stop("EM likelihood decreased")
      history <- c(history,next_fit$loglik)
      current <- next_fit; intercept <- new_c; phi <- new_phi
      sigma <- new_sigma; P <- new_P
      if(abs(gain)<tol*(1+abs(current$loglik))) {converged<-TRUE;break}
    }
    candidates[[start]] <- list(fit=current,history=history,converged=converged,
                                iterations=length(history)-1L)
  }
  valid <- which(vapply(candidates,function(x)x$converged,logical(1)))
  if(!length(valid)) stop("No EM start converged; increase maxit or revise specification")
  chosen <- valid[which.max(vapply(candidates[valid],function(x)x$fit$loglik,numeric(1)))]
  best <- candidates[[chosen]]
  # Labels: low intercept first, or low variance first for the zero-mean model.
  order <- if(zero_mean) order(best$fit$sigma) else order(best$fit$c)
  f <- best$fit
  f <- ms_filter(obs,f$c[order],f$phi[order],f$sigma[order],f$P[order,order],
                 y0,prior[order])
  f$index <- index; f$converged <- TRUE; f$history <- best$history
  f$iterations <- best$iterations; f$successful_starts <- length(valid)
  f$start_loglik <- vapply(candidates,function(x)x$fit$loglik,numeric(1))
  f$sd_floor <- sd_floor; f$probability_floor <- probability_floor
  f$boundary <- min(f$sigma)<=sd_floor*(1+1e-6) ||
    any(f$P<=probability_floor*(1+1e-6)) ||
    (ar && max(abs(f$phi))>=phi_bound*(1-1e-6))
  f$npar <- (if(zero_mean) 0L else 2L)+as.integer(ar)+
    (if(switching_variance) 2L else 1L)+2L
  f$aic <- -2*f$loglik+2*f$npar
  f$bic <- -2*f$loglik+log(n)*f$npar
  f
}

# Exact conditional first and second moments, integrating over future regimes.
#' Exact moments of a regime-uncertain AR(1) forecast
#'
#' Propagate joint state probabilities and state-weighted first and second moments to integrate over all future regimes and innovations.
#'
#' @param c Numeric vector of regime-specific recursion intercepts, not unconditional regime means.
#' @param phi Numeric vector of regime-specific AR(1) coefficients.
#' @param sigma Positive numeric vector of regime innovation standard deviations.
#' @param P Square row-stochastic transition matrix: P[i,j] is the probability of moving from state i to state j. Probability vectors are columns.
#' @param last_y Finite scalar observed at the forecast origin.
#' @param state_prob Numeric vector of filtered state probabilities at the origin.
#' @param horizon Positive integer number of future observations.
#' @return A list containing horizon-length mean and variance vectors for future observations, and a horizon-by-K state_probability matrix.
#' @details Parameters are fixed and transitions are autonomous. At horizons beyond one, the future lagged observation and regime need not be conditionally independent. The calculation retains their joint moments instead of replacing them by products of marginal means. Forecast variance includes both innovation and regime uncertainty; it excludes estimation uncertainty. A Gaussian prediction interval need not be valid for the resulting mixture.
#' @examples
#' P <- matrix(c(.9,.1,.05,.95),2,byrow=TRUE)
#' f <- ms_forecast_moments(0,c(.4,.6),c(-1,1),c(.3,.3),c(.5,.5),P,horizon=6)
#' cbind(mean=f$mean,sd=sqrt(f$variance))
#' stopifnot(all(f$variance > 0))
#' @family Markov-switching models
#' @export
ms_forecast_moments <- function(last_y, state_prob, c, phi, sigma, P, horizon=12L) {
  K <- ms_check_model(c,phi,sigma,P)
  if(!.idss_scalar(last_y) || !.idss_prob(state_prob,K) || !.idss_integer(horizon)) stop("Invalid forecast inputs")
  a <- state_prob; b <- last_y*a; d <- last_y^2*a
  means <- variances <- numeric(horizon); states <- matrix(0,horizon,K)
  for(h in seq_len(horizon)) {
    next_a <- as.vector(crossprod(P,a)); next_b <- as.vector(crossprod(P,b))
    next_d <- as.vector(crossprod(P,d))
    b <- c*next_a+phi*next_b
    d <- (c^2+sigma^2)*next_a+2*c*phi*next_b+phi^2*next_d
    a <- next_a
    means[h] <- sum(b); variances[h] <- sum(d)-sum(b)^2; states[h,] <- a
  }
  list(mean=means,variance=variances,state_probability=states)
}

# Two-state simulation of the predictive distribution; fitted parameters fixed.
#' Simulate a two-state switching-AR predictive distribution
#'
#' Generate future observation paths conditional on a filtered two-state distribution and fixed model parameters.
#'
#' @param c Numeric vector of regime-specific recursion intercepts, not unconditional regime means.
#' @param phi Numeric vector of regime-specific AR(1) coefficients.
#' @param sigma Positive numeric vector of regime innovation standard deviations.
#' @param P Square row-stochastic transition matrix: P[i,j] is the probability of moving from state i to state j. Probability vectors are columns.
#' @param last_y Finite scalar observed at the forecast origin.
#' @param state_prob Numeric vector of filtered state probabilities at the origin.
#' @param horizon Positive integer number of future observations.
#' @param draws Positive integer number of simulated paths.
#' @return Numeric matrix with draws rows and horizon columns; each row is one future path.
#' @details This simulator supports exactly two regimes. It first draws the current state from state_prob, transitions using P, and draws independent future Gaussian innovations. Quantiles at each horizon give pointwise mixture prediction intervals, not simultaneous path bands. Parameter-estimation uncertainty is not simulated.
#' @examples
#' set.seed(123)
#' P <- matrix(c(.9,.1,.05,.95),2,byrow=TRUE)
#' x <- ms_forecast_draws(0,c(.4,.6),c(-1,1),c(.3,.3),c(.5,.5),P,
#'                        horizon=4,draws=2000)
#' apply(x,2,quantile,probs=c(.05,.5,.95))
#' @family Markov-switching models
#' @export
ms_forecast_draws <- function(last_y, state_prob, c, phi, sigma, P,
                              horizon=12L, draws=20000L) {
  if(ms_check_model(c,phi,sigma,P)!=2L) stop("This simulator uses two states")
  if(!.idss_scalar(last_y) || !.idss_prob(state_prob,2L) ||
     !.idss_integer(horizon) || !.idss_integer(draws))
    stop("Invalid forecast inputs")
  state <- ifelse(runif(draws)<state_prob[1],1L,2L)
  y <- rep(last_y,draws); result <- matrix(0,draws,horizon)
  for(h in seq_len(horizon)) {
    state <- ifelse(runif(draws)<P[state,1],1L,2L)
    y <- c[state]+phi[state]*y+sigma[state]*rnorm(draws)
    result[,h] <- y
  }
  result
}
