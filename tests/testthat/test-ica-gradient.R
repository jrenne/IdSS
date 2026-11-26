test_that("ICA gradients agree with finite differences away from the identity", {
  set.seed(20261126)
  for (n in 2:4) {
    errors <- matrix(stats::rt(240*n,df=7),240,n)
    Y <- scale(errors,center=TRUE,scale=FALSE) %*% solve(chol(stats::cov(errors)))
    count <- n*(n-1)/2
    for (type in c("student","gaussian","mixt.gaussian","hyper.sec","laplace")) {
      distri <- list(type=rep(type,n),df=rep(7,n),p=rep(.4,n),
                     mu=rep(.2,n),sigma=rep(.7,n))
      for (theta in list(rep(0,count),seq_len(count)*.07-.12)) {
        actual <- d.func.2.minimize(theta,Y,distri)
        h <- 1e-6
        numerical <- vapply(seq_len(count),function(i) {
          delta <- rep(0,count);delta[i] <- h
          (func.2.minimize(theta+delta,Y,distri)-
           func.2.minimize(theta-delta,Y,distri))/(2*h)
        },numeric(1))
        expect_length(actual,count)
        expect_true(all(is.finite(actual)))
        expect_lt(max(abs(actual-numerical)),1e-5,
                  label=paste("gradient discrepancy",n,type,paste(theta,collapse=",")))
      }
    }
  }
})
