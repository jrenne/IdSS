# Shared internal validation for latent-state procedures.
.idss_real <- function(x) is.numeric(x) && !is.complex(x)
.idss_scalar <- function(x) .idss_real(x) && length(x)==1L && is.finite(x)
.idss_integer <- function(x,lower=1L) .idss_scalar(x) && x>=lower &&
  x<=.Machine$integer.max && x==floor(x)
.idss_flag <- function(x) is.logical(x) && length(x)==1L && !is.na(x)
.idss_prob <- function(x,k) .idss_real(x) && is.null(dim(x)) && length(x)==k &&
  all(is.finite(x)) && all(x>=0) && abs(sum(x)-1)<=1e-10
