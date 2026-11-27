# Independent mathematical checks motivated by the companion-course review.
sample_var <- function(T=220,n=3) {
  Y <- matrix(stats::rnorm(T*n),T,n)
  for(t in 2:T) Y[t,] <- .35*Y[t-1,]+Y[t,]
  colnames(Y) <- paste0("y",seq_len(n));Y
}

test_that("LP impact responses match Cholesky and do not depend on maximum horizon", {
  set.seed(2801);Y <- sample_var()
  fit <- vars::VAR(Y,p=1);B <- t(chol(stats::cov(stats::residuals(fit))))
  for(k in 1:3) {
    short <- suppressMessages(make.jorda.irf(Y,k,1,1,indic.plot=0))
    long <- suppressMessages(make.jorda.irf(Y,k,5,1,indic.plot=0))
    expect_equal(as.numeric(short$IRFs[1,]),B[,k],tolerance=1e-9,ignore_attr=TRUE)
    expect_equal(short$IRFs,long$IRFs[1:2,],tolerance=1e-12)
  }
  set.seed(2801);a <- suppressMessages(make.jorda.irf(Y,2,4,1,indic.plot=0))
  Y[,2] <- 7*Y[,2]
  b <- suppressMessages(make.jorda.irf(Y,2,4,1,indic.plot=0))
  expect_equal(b$IRFs[,c(1,3)],a$IRFs[,c(1,3)],tolerance=1e-9)
  expect_equal(b$IRFs[,2],7*a$IRFs[,2],tolerance=1e-9)
})

test_that("LP HAC covariance uses the actual horizon sample size", {
  set.seed(2802);Y <- sample_var(120,2);h <- 4;k <- 2;p <- 1;q <- 2
  ans <- suppressMessages(make.jorda.irf(Y,k,h,p,indic.plot=0,hac.lag=q))
  dates <- (p+1):(nrow(Y)-h);X <- cbind(1,Y[dates,],Y[dates-1,])
  b <- qr.solve(X,Y[dates+h,1]);scores <- X*as.numeric(Y[dates+h,1]-X%*%b)
  scores <- scale(scores,center=TRUE,scale=FALSE)
  meat <- crossprod(scores)
  for(j in 1:q) {
    lagged <- crossprod(scores[(j+1):nrow(scores),,drop=FALSE],scores[1:(nrow(scores)-j),,drop=FALSE])
    meat <- meat+(1-j/(q+1))*(lagged+t(lagged))
  }
  inverse <- solve(crossprod(X));V <- inverse%*%meat%*%inverse
  full <- cbind(1,Y[2:nrow(Y),1],Y[1:(nrow(Y)-1),])
  shock.sd <- stats::sd(stats::lm.fit(full,Y[2:nrow(Y),2])$residuals)
  expect_equal(ans$IRFs[h+1,1],as.numeric(b[3]*shock.sd),tolerance=1e-10)
  expect_equal(ans$all.stdv.IRFs[h+1,1],sqrt(V[3,3])*shock.sd,tolerance=1e-10)
})

test_that("LP-IV has unit impact and accepts explicit bandwidths", {
  set.seed(2803);z <- stats::rnorm(180);u <- stats::rnorm(180)
  y <- z+u;Y <- cbind(policy=y,outcome=.7*y+.5*u+stats::rnorm(180))
  a <- suppressMessages(make.LPIV.irf(Y,matrix(z,ncol=1),3,1,1,0,hac.lag=1))
  b <- suppressMessages(make.LPIV.irf(Y,matrix(z,ncol=1),3,1,1,0,hac.lag=4))
  expect_equal(a$IRFs[1,1],1,tolerance=1e-10)
  expect_equal(a$IRFs,b$IRFs,tolerance=1e-12)
  expect_true(all(is.finite(a$all.stdv.IRFs)))
  expect_gt(max(abs(a$all.stdv.IRFs-b$all.stdv.IRFs)),1e-5)
})

test_that("FEV maximization windows retain early response contributions", {
  set.seed(2804);Y <- sample_var(260,3)
  baseline <- IdSS:::svar.ordering.all(Y,2,8,3)
  S <- matrix(0,3,3)
  for(H in 3:6) for(j in seq_len(H)) S <- S+tcrossprod(baseline$IRFs[1,,j])
  ans <- svar.fevmax(Y,2,1,"test shock",3,6,1,1,8,1,indic.plot=0)
  Q <- solve(baseline$B.hat,ans$simulated.IRFs[,,1,1])
  expect_equal(as.numeric(crossprod(Q[,1],S%*%Q[,1])),max(eigen(S,symmetric=TRUE)$values),tolerance=1e-9)
  expect_equal(crossprod(Q),diag(3),tolerance=1e-10)
})

test_that("Monte Carlo completes every genuine covariance draw after rejection", {
  set.seed(2805);Y <- sample_var(9,2)
  baseline <- svar.ordering(Y,1,1,3,inference=0,indic.plot=0)$IRFs
  set.seed(2806);ans <- svar.ordering(Y,1,1,3,inference=3,nb.draws=50,indic.plot=0)
  expect_equal(dim(ans$simulated.IRFs),c(3L,2L,50L))
  expect_true(all(is.finite(ans$simulated.IRFs)))
  equals.baseline <- vapply(1:50,function(i) isTRUE(all.equal(ans$simulated.IRFs[,,i],baseline,tolerance=1e-13)),logical(1))
  expect_false(any(equals.baseline))
})

sign_inputs <- function(n) list(signs=rep(list(matrix(0,n,n)),n),h=rep(list(1),n))

test_that("QR and recursive unrestricted rotations have centered Haar margins", {
  set.seed(2807);Y <- sample_var(100,3);s <- sign_inputs(3)
  C <- IdSS:::svar.ordering.all(Y,1,2,3)$B.hat
  for(method in 0:1) {
    ans <- svar.signs(Y,1,3,nb.periods.IRF=2,nb.draws=1000,indic.plot=0,
                     sign.restrictions=s$signs,horizon=s$h,recursive=method)
    expect_equal(ans$xx,1000)
    rotations <- vapply(1:1000,function(i) as.numeric(solve(C,ans$IRFs.signs[,,1,i])),numeric(9))
    expect_lt(max(abs(rowMeans(rotations))),.09)
    expect_lt(max(abs(rowMeans(rotations^2)-1/3)),.07)
    dets <- apply(rotations,2,function(q) det(matrix(q,3,3)))
    expect_gt(mean(dets>0),.4);expect_lt(mean(dets>0),.6)
    for(j in 1:3) {
      expected <- t(apply(ans$IRFs.signs[,j,,,drop=FALSE],c(1,3),stats::median))
      expect_equal(unname(ans$all.CI.median[[j]]),unname(expected),tolerance=1e-12)
    }
  }
})

test_that("every parameter bootstrap enters sign-restricted draws", {
  set.seed(2808);Y <- sample_var(100,2);s <- sign_inputs(2)
  ans <- svar.signs(Y,1,2,nb.periods.IRF=3,bootstrap.replications=3,
                   nb.draws=7,indic.plot=0,sign.restrictions=s$signs,horizon=s$h)
  expect_equal(ans$xx,21)
  for(i in 1:21) {
    group <- (i-1)%/%7+1
    C <- ans$simulated.IRFs[,,1,group]
    Q <- solve(C,ans$IRFs.signs[,,1,i])
    expect_equal(crossprod(Q),diag(2),tolerance=1e-9)
  }
})

test_that("zero restrictions use rank and contradictory inequalities are rejected", {
  set.seed(2809);Y <- sample_var(100,3);s <- sign_inputs(3)
  zero <- list(rbind(c(1,0,0),c(1,0,0)))
  ans <- svar.signs(Y,1,1,nb.periods.IRF=2,nb.draws=40,indic.plot=0,
                   sign.restrictions=s$signs,horizon=s$h,SR.restrictions=zero)
  expect_equal(ans$xx,40)
  expect_lt(max(abs(ans$IRFs.signs[1,1,1,])),1e-10)
  s$signs[[1]] <- rbind(c(1,0,0),c(-1,0,0))
  expect_error(svar.signs(Y,1,1,nb.periods.IRF=2,nb.draws=30,indic.plot=0,
                        sign.restrictions=s$signs,horizon=s$h),"No admissible rotations")
})

test_that("bias-corrected bootstrap is invariant to translations of the data", {
  set.seed(2810);Y <- sample_var(100,2)
  set.seed(2811);a <- svar.ordering(Y,1,1,5,inference=4,nb.draws=16,indic.plot=0)
  shifted <- sweep(Y,2,c(100,-40),"+")
  set.seed(2811);b <- svar.ordering(shifted,1,1,5,inference=4,nb.draws=16,indic.plot=0)
  expect_equal(a$IRFs,b$IRFs,tolerance=1e-8)
  expect_equal(a$simulated.IRFs,b$simulated.IRFs,tolerance=1e-7)
  expect_true(all(is.finite(a$simulated.IRFs)))
})
