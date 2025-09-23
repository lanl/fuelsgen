#' @title Plot MCMC posteriors
#'
#' @description Pairs plot of posterior density estimates and prior distributions
#' @param theta.list list of matrices containing posterior samples
#' @param truth optional vector of known parameter values
#' @param cols colors for posterior density estimate plots
#' @param labels legend names for posterior density estimates
#' @export
#'
contour_pairs <- function(theta.list, prior, truth = NULL,
                          cols = c('cornflowerblue','darkorange','forestgreen'),
                          rho.bw = "nrd0", mu.bw = "nrd0", var.bw = "nrd0",
                          lambda.bw = "nrd0", nu.bw = "nrd0",
                          labels = c('posterior')) {
  
  stopifnot(length(theta.list) >= 1)
  
  # ---- pull prior component objects (used for bounds AND prior curves)
  rho.prior <- prior$prior_params[[which(prior$names == 'rho')]]
  mu.prior  <- prior$prior_params[[which(prior$names == 'mu')]]
  var.prior <- prior$prior_params[[which(prior$names == 's2')]]
  lam.prior <- prior$prior_params[[which(prior$names == 'lambda')]]
  nu.prior  <- prior$prior_params[[which(prior$names == 'nu')]]
  
  # ---- build prior curves on their own grids (store back into objects)
  rho.prior$x <- seq(rho.prior$bounds[1], rho.prior$bounds[2], length.out = 1000)
  rho.prior$y <- dunif(rho.prior$x, min = rho.prior$bounds[1], max = rho.prior$bounds[2])
  
  mu.prior$x  <- seq(mu.prior$bounds[1], mu.prior$bounds[2], length.out = 1000)
  mu.prior$y  <- truncdist::dtrunc(mu.prior$x, 'norm',
                                   a = mu.prior$params[3], b = mu.prior$params[4],
                                   mean = mu.prior$params[1], sd = mu.prior$params[2])
  
  var.prior$x <- seq(var.prior$bounds[1], var.prior$bounds[2], length.out = 1000)
  var.prior$y <- LaplacesDemon::dhalfcauchy(var.prior$x, scale = var.prior$params[1])
  
  lam.prior$x <- seq(lam.prior$bounds[1], lam.prior$bounds[2], length.out = 1000)
  lam.prior$y <- dgamma(lam.prior$x, shape = lam.prior$params[1], scale = lam.prior$params[2])
  
  nu.prior$x  <- seq(nu.prior$bounds[1], nu.prior$bounds[2], length.out = 1000)
  nu.prior$y  <- dgamma(nu.prior$x, shape = nu.prior$params[1], scale = nu.prior$params[2])
  
  # ---- mapping: columns in theta to parameter names & transforms
  # theta cols: 1=rho, 2=mu, 3=sigma (plot as sigma^2), 4=lambda, 5=nu
  par_names  <- c("rho","mu","s2","lambda","nu")
  priors     <- list(rho.prior, mu.prior, var.prior, lam.prior, nu.prior)
  bw_list    <- list(rho.bw,    mu.bw,    var.bw,    lambda.bw,  nu.bw)
  
  get_col <- function(mat, j) {
    mat[,j]
  }
  
  # how many posterior runs to overlay on each subplot
  n.post <- length(theta.list)
  n.par  <- length(par_names)  # = 5
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(n.par, n.par), mar = c(2,2,.5,.5))
  
  # helper: draw prior+posteriors marginal for parameter j (diagonal)
  draw_marginal <- function(j) {
    pr  <- priors[[j]]
    bw  <- bw_list[[j]]
    
    # posterior densities per run, clipped to bounds
    post_dens <- vector("list", n.post)
    ymax <- 0
    for (i in 1:n.post) {
      xj <- get_col(theta.list[[i]], j)
      post_dens[[i]] <- density(xj,
                                from = pr$bounds[1],
                                to   = pr$bounds[2],
                                bw   = bw)
      ymax <- max(ymax, max(post_dens[[i]]$y, na.rm=TRUE))
    }
    
    # plot prior
    plot(pr$x, pr$y, type = 'l', lty = 2,
         xlim = pr$bounds, ylim = c(0, max(ymax, max(pr$y, na.rm=TRUE))),
         xlab = '', ylab = '', yaxt='n')
    # overlay posteriors
    for (i in 1:n.post) lines(post_dens[[i]], col = cols[i])
    
    # optional truth line
    if (!is.null(truth) && length(truth) >= j) {
      v <- if (j == 3) truth[3] else truth[j]   # truth[3] already variance in your original
      abline(v = v, col = 'red', lty = 2)
    }
    
    # tiny corner label
    lbl <- switch(par_names[j],
                  rho = expression(rho),
                  mu  = expression(mu),
                  s2  = expression(sigma^2),
                  lambda = expression(lambda),
                  nu = expression(nu))
    legend(inset=c(.05,.05),'topright', legend = lbl, bty='n', cex=1.2)
  }
  
  # helper: 2D KDE contours between params (r != c)
  draw_bivar <- function(r, c) {
    pr_r <- priors[[r]]
    pr_c <- priors[[c]]
    first <- TRUE
    for (i in 1:n.post) {
      xr <- get_col(theta.list[[i]], r)
      xc <- get_col(theta.list[[i]], c)
      kd <- MASS::kde2d(xr, xc)
      if (first) {
        contour(kd, nlevels = 10, col = cols[i],
                xlim = pr_r$bounds, ylim = pr_c$bounds,
                xlab = '', ylab = '')
        first <- FALSE
      } else {
        contour(kd, nlevels = 10, col = cols[i], add = TRUE)
      }
    }
    if (!is.null(truth) && length(truth) >= max(r,c)) {
      tr <- if (r == 3) truth[3] else truth[r]
      tc <- if (c == 3) truth[3] else truth[c]
      points(tr, tc, col = 'red', pch = 16)
    }
  }
  
  # ---- draw the full n.par x n.par grid
  for (r in 1:n.par) {
    for (c in 1:n.par) {
      if (r == c) {
        draw_marginal(r)
      } else if(r>c){
        draw_bivar(c,r)
      } else if(c==2 & r==1){
        # legend panel: top-left empty cell overlayed with legend (optional)
        # par(xpd = NA)  # allow drawing over margins
        plot.new()
        if (!is.null(labels) && length(labels) >= n.post) {
          # create a small legend box; this will appear after the grid in many devices
          legend('center',
                 legend = c('prior', labels[seq_len(n.post)], if (!is.null(truth)) 'truth' else NULL),
                 lty    = c(2, rep(1, n.post), if (!is.null(truth)) 2 else NULL),
                 col    = c('black', cols[seq_len(n.post)], if (!is.null(truth)) 'red' else NULL),
                 bty='n', cex=1.1)
        }
      } else{
        plot.new()
      }
    }
  }
}

#' @title Plot samples of MCMC posteriors
#'
#' @description Pairs plot of posterior density estimates and prior distributions
#' @param theta.list list of matrices containing posterior samples
#' @param truth optional vector of known parameter values
#' @param cols colors for posterior density estimate plots
#' @param labels legend names for posterior density estimates
#' @export
#'
scatter_pairs <- function(theta.list, prior, truth = NULL,
                          cols = c('cornflowerblue','darkorange','forestgreen'),
                          rho.bw = "nrd0", mu.bw = "nrd0", var.bw = "nrd0",
                          lambda.bw = "nrd0", nu.bw = "nrd0",
                          labels = c('posterior'),
                          pt.max = 3000,       # max points per chain per panel
                          pt.alpha = 0.25,     # point transparency
                          pt.cex = 0.4,        # point size
                          pt.pch = 16,         # point shape
                          lwd = 2) {       
  stopifnot(length(theta.list) >= 1)
  
  # ---- pull prior component objects
  rho.prior <- prior$prior_params[[which(prior$names == 'rho')]]
  mu.prior  <- prior$prior_params[[which(prior$names == 'mu')]]
  var.prior <- prior$prior_params[[which(prior$names == 's2')]]
  lam.prior <- prior$prior_params[[which(prior$names == 'lambda')]]
  nu.prior  <- prior$prior_params[[which(prior$names == 'nu')]]
  
  # ---- prior curves on their grids
  rho.prior$x <- seq(rho.prior$bounds[1], rho.prior$bounds[2], length.out = 1000)
  rho.prior$y <- dunif(rho.prior$x, min = rho.prior$bounds[1], max = rho.prior$bounds[2])
  
  mu.prior$x  <- seq(mu.prior$bounds[1], mu.prior$bounds[2], length.out = 1000)
  mu.prior$y  <- truncdist::dtrunc(mu.prior$x, 'norm',
                                   a = mu.prior$params[3], b = mu.prior$params[4],
                                   mean = mu.prior$params[1], sd = mu.prior$params[2])
  
  var.prior$x <- seq(var.prior$bounds[1], var.prior$bounds[2], length.out = 1000)
  var.prior$y <- LaplacesDemon::dhalfcauchy(var.prior$x, scale = var.prior$params[1])
  
  lam.prior$x <- seq(lam.prior$bounds[1], lam.prior$bounds[2], length.out = 1000)
  lam.prior$y <- dgamma(lam.prior$x, shape = lam.prior$params[1], scale = lam.prior$params[2])
  
  nu.prior$x  <- seq(nu.prior$bounds[1], nu.prior$bounds[2], length.out = 1000)
  nu.prior$y  <- dgamma(nu.prior$x, shape = nu.prior$params[1], scale = nu.prior$params[2])
  
  # ---- mapping / helpers
  par_names <- c("rho","mu","s2","lambda","nu")
  priors    <- list(rho.prior, mu.prior, var.prior, lam.prior, nu.prior)
  bw_list   <- list(rho.bw,    mu.bw,    var.bw,    lambda.bw,  nu.bw)
  
  get_col <- function(mat, j) mat[, j]
  
  n.post <- length(theta.list)
  n.par  <- length(par_names)
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(n.par, n.par), mar = c(2,2,.5,.5))
  
  # diagonal: marginal prior + posterior densities
  draw_marginal <- function(j) {
    pr <- priors[[j]]; bw <- bw_list[[j]]
    post_dens <- vector("list", n.post)
    ymax <- 0
    for (i in 1:n.post) {
      xj <- get_col(theta.list[[i]], j)
      if(diff(range(xj))==0){
        post_dens[[i]] = list(x=xj[1],y=0)
      } else{
        post_dens[[i]] <- density(xj,
                                  from = pr$bounds[1],
                                  to   = pr$bounds[2],
                                  bw   = bw)
      }
      ymax <- max(ymax, max(post_dens[[i]]$y, na.rm=TRUE))
    }
    plot(pr$x, pr$y, type='l', lty=2,
         xlim = pr$bounds, ylim = c(0, max(ymax, max(pr$y, na.rm=TRUE))),
         xlab='', ylab='', yaxt='n')
    for (i in 1:n.post){
      if(length(post_dens[[i]]$x)==1){ # fixed parameter
        abline(v = post_dens[[i]]$x, col = cols[i], lwd=lwd)
      } else{
        lines(post_dens[[i]], col = cols[i],lwd=lwd)  
      }
    }
    if (!is.null(truth) && length(truth) >= j) {
      v <- truth[j]
      abline(v = v, col = 'red', lty = 2)
    }
    lbl <- switch(par_names[j],
                  rho = expression(rho),
                  mu  = expression(mu),
                  s2  = expression(sigma^2),
                  lambda = expression(lambda),
                  nu = expression(nu))
    legend(inset=c(.05,.05),'topright', legend=lbl, bty='n', cex=1.2)
  }
  
  # lower triangle: scatter instead of contours
  draw_bivar_scatter <- function(r, c) {
    pr_r <- priors[[r]]; pr_c <- priors[[c]]
    # set up empty plot with bounds
    plot(NA, NA, xlim = pr_r$bounds, ylim = pr_c$bounds,
         xlab = '', ylab = '', xaxt='s', yaxt='s')
    for (i in 1:n.post) {
      xr <- get_col(theta.list[[i]], r)
      xc <- get_col(theta.list[[i]], c)
      n <- length(xr)
      take <- if (n > pt.max) sample.int(n, pt.max) else seq_len(n)
      points(xr[take], xc[take],
             col = grDevices::adjustcolor(cols[i], alpha.f = pt.alpha),
             pch = pt.pch, cex = pt.cex)
    }
    if (!is.null(truth) && length(truth) >= max(r,c)) {
      points(truth[r], truth[c], col='red', pch=16)
    }
  }
  
  # draw full grid: diag marginals; lower scatter; (1,2) legend; others blank
  for (r in 1:n.par) {
    for (c in 1:n.par) {
      if (r == c) {
        draw_marginal(r)
      } else if (r > c) {
        # keep axis orientation consistent with cell (r,c)
        draw_bivar_scatter(c, r)
      } else if (c == 2 && r == 1) {
        plot.new()
        if (!is.null(labels) && length(labels) >= n.post) {
          legend('center',
                 legend = c('prior', labels[seq_len(n.post)], if (!is.null(truth)) 'truth' else NULL),
                 lty    = c(2, rep(1, n.post), if (!is.null(truth)) 2 else NULL),
                 col    = c('black', cols[seq_len(n.post)], if (!is.null(truth)) 'red' else NULL),
                 bty='n', cex=1.1)
        }
      } else {
        plot.new()
      }
    }
  }
}
