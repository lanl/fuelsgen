# --- log-weights for CMP pmf ---
.cmp_logw <- function(k, lambda, nu) {
  k * log(lambda) - nu * lgamma(k + 1)
}
# --- Build pmf on 0..K with a rigorous relative tail bound and stable normalization ---
cmp_pmf_trunc_strict <- function(lambda, nu, rel_tail = 1e-14,
                                 K_init = NULL, K_cap = 1e6) {
  stopifnot(is.finite(lambda), is.finite(nu), lambda > 0, nu > 0)
  # start near the mode; floor(lambda^(1/nu)) is a good anchor
  K <- if (is.null(K_init)) max(50L, floor(lambda^(1/nu)) + 10L) else as.integer(K_init)
  
  repeat {
    k  <- 0:K
    lw <- .cmp_logw(k, lambda, nu)
    m  <- max(lw)
    v  <- exp(lw - m)         # unnormalized weights (scaled)
    S  <- sum(v)              # sum of scaled weights
    # ratio of successive terms beyond K
    r_next <- lambda / ((K + 1)^nu)
    
    # only assess the geometric tail bound once we're past the mode (r_next < 1)
    if (r_next < 1) {
      # relative tail mass bound (w_K/S) * r_next/(1 - r_next)
      vK <- v[length(v)]
      tail_rel <- (vK / S) * (r_next / (1 - r_next))
      if (is.finite(tail_rel) && tail_rel < rel_tail) {
        p   <- v / S
        cdf <- cumsum(p)
        return(list(p = p, cdf = cdf, K = K, mode = which.max(lw) - 1L))
      }
    }
    
    # otherwise, extend support
    K <- as.integer(min(K_cap, ceiling(1.25 * K) + 20L))
    if (K >= K_cap) stop("cmp_pmf_trunc_strict: K_cap reached; increase cap or check parameters.")
  }
}

# --- Vectorized inverse-CDF sampler (no Poisson shortcut at nu=1) ---
cmp_sampler_factory <- function(lambda, nu, rel_tail = 1e-14) {
  pmf <- cmp_pmf_trunc_strict(lambda, nu, rel_tail = rel_tail)
  cdf <- pmf$cdf
  function(n) {
    u <- runif(n)
    as.integer(findInterval(u, cdf))
  }
}

# --- Newton solve for lambda given (mu, nu), using strict pmf moments ---
mu_nu_to_lambda_newton <- function(mu, nu, tol = 1e-10, maxit = 50,
                                   rel_tail = 1e-14) {
  stopifnot(is.finite(mu), is.finite(nu), mu >= 0, nu > 0)
  # generic initial guess; for nu≈1 this is already close
  lam <- max(1e-12, mu^nu)
  
  for (it in 1:maxit) {
    mm   <- cmp_moments_from_pmf(lam, nu, rel_tail = rel_tail)
    EN   <- mm$mean
    VarN <- mm$var
    if (!is.finite(EN) || !is.finite(VarN)) stop("Non-finite moment in Newton step.")
    if (abs(EN - mu) <= tol * max(1, mu)) return(lam)
    
    # dE/dlambda = Var/λ (score identity under exponential family weights)
    grad <- VarN / max(lam, 1e-12)
    step <- (EN - mu) / max(grad, 1e-16)
    
    lam_new <- lam - step
    if (!is.finite(lam_new) || lam_new <= 0) {
      lam <- lam * 0.5   # backtrack if needed
    } else {
      lam <- lam_new
    }
  }
  warning("mu_nu_to_lambda_newton: max iterations reached; returning last lambda.")
  lam
}

# --- Moments from the same (strict) pmf ---
cmp_moments_from_pmf <- function(lambda, nu, rel_tail = 1e-14) {
  pmf <- cmp_pmf_trunc_strict(lambda, nu, rel_tail = rel_tail)
  k <- 0:pmf$K
  p <- pmf$p
  EN  <- sum(k * p)
  EN2 <- sum(k * k * p)
  list(mean = EN, var = EN2 - EN * EN)
}

# --- Public: strict mean-parameterized CMP sampler ---
my_rcomp <- function(n, mu, nu, rel_tail = 1e-14) {
  if (length(n) > 1) n <- length(n)
  stopifnot(length(mu) == 1L, length(nu) == 1L)
  if (mu == 0) return(integer(n))
  
  lambda <- mu_nu_to_lambda_newton(mu, nu, rel_tail = rel_tail)
  sampler <- cmp_sampler_factory(lambda, nu, rel_tail = rel_tail)
  sampler(n)
}