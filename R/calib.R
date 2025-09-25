#' @title MCMC Calibration for FuelGen model
#'
#' @description Calibrate parameters to observed data using addaptive metropolis
#' @param y_obs metrics for observed data
#' @param fuelsgen object for observed data
#' @param prior prior from get_prior_info()
#' @param adapt.par adaptive metropolis parameters. The first is the number of iterations to do before starting adaptation. The second is how often to adapt. The third is the proportion of previous samples to use for adaptation. The fourth is the proportion of samples after which to stop adaptation.
#' @param prop.sigma initial proposal covariance for adaptive metropolis
#' @param GP.init.size grid size for discretization of intensity function. Larger gives a more accurate approximation, but increases computational complexity
#' @param I.transform 'exp' or 'logistic', the transormation used to convert the GP realizations to the positive real line
#' @param gen_reps the number of simulated realizations for each proposed value of the model parameters. More reps gives more robust sampling algorithm, but can be computationally prohibitive.
#' @param ABC use ABC sampling rather than LLH based inference
#' @param ABC_eps distance cutoff used for ABC, needs to be tuned based on the metrics used
#' @param ABC_wt weight used in ABC distance metric. The second element is for the summary statistics defining the mean and variance of the number of shrubs over replicate observations, the first is for all the other metrics.
#' @param gen_parallel generate realizations in parallel - often slower than in series
#' @param mets_parallel compute metrics in parallel - usually faster than in series
#' @param make_cluster make the parallel computing cluster within the function call
#' @param verbose print status
#' @param fixed named vector of fixed parameters (ex. c('nu'=1,'mu'=1.5) fixes dispersion at 1 and radius mean at 1.5)
#' @export
#'
mcmc_MH_adaptive = function(y_obs,fuel,prior,
                            n.samples=10000,n.burn=1000,
                            adapt.par = c(100,20,.5,.75),
                            prop.sigma = diag(.1^2*prior$theta_est), 
                            GP.init.size = 32, I.transform = 'exp',
                            gen_reps = 25, 
                            ABC = F, ABC_eps = 10, ABC_wt = c(1,.1),
                            gen_parallel = F, mets_parallel = T,
                            make_cluster=T,verbose=T,fixed=NULL)
{
  llh_mh_adaptive = function(theta, ABC, ABC_eps, ABC_wt)
  {
    # if we specified that any parameters should be fixed, do it here
    if(!is.null(fixed)){
      for(i in 1:length(fixed)){
        fixed_idx <- match(names(fixed)[i], prior$names)
        theta[fixed_idx] = fixed[i]
      }
    }
    if(any(theta<prior$lb) | any(theta>prior$ub)){
      if(ABC){
        return(Inf)
      } else{
        return(-Inf)
      }
    }
    sim_fuels = gen_data(theta,fuel$dimX,fuel$dimY,1,fuel$X.locs,fuel$X.vals,fuel$Beta,gen_reps,GP.init.size,NULL,I.transform,.217622,F)
    
    metrics = get_mets(sim_fuels, y_obs$info, mets_parallel, make_cluster = F)
    metrics = fix_nan_metrics(theta, metrics, fuel)
  
    if(ABC){
      D = abc_distance(y_obs$mets, metrics$mets, w_met = ABC_wt[1], w_cnt = ABC_wt[2])
      return(D)
    } else{
      ##################
      # precondition the metrics using the eigen decomp of the prior metrics covariance
      # seems to work really well
      
      d <- ncol(prior$Sigma)
      
      # Eigen "surgery": floor tiny eigenvalues before whitening
      e <- eigen((prior$Sigma + t(prior$Sigma))/2, symmetric=TRUE)
      lam  <- pmax(e$values, median(e$values) * 1e-3)  # floor at 1e-3 * median
      Linv <- t(e$vectors %*% diag(1/sqrt(lam)) %*% t(e$vectors))  # whitening matrix
      
      # Whiten both observed and simulated metrics BEFORE likelihood
      Yw <- y_obs$mets %*% Linv
      
      # In whitened space, use prior Sigma = I_d with moderate nu
      prior$Sigma <- diag(d)
      prior$nu    <- max(prior$nu, d + 5)
      
      Mw <- metrics$mets %*% Linv
      Sigma_gen = cov(Mw)
      Psi_n = lqmm::make.positive.definite((prior$nu-1)*prior$Sigma + (gen_reps-1)*Sigma_gen)
      Sigma = Psi_n / (prior$nu + gen_reps + ncol(Psi_n) + 1) # mode of invWishart
      ##################
      
      Sinv = chol2inv(chol(Sigma))
      ldetS = determinant(Sigma)$modulus
      
      # average over the metrics and compute the likelihood once
      metrics = colMeans(metrics$mets)
      llh = numeric(nrow(y_obs$mets))
      for(i in 1:length(llh)){
        # loop over replicate observations
        llh[i] = - 0.5 * (ldetS + ((y_obs$mets[i,,drop=F] - metrics) %*% Sinv %*% t(y_obs$mets[i,,drop=F] - metrics)))
      }
      return(sum(llh)+lprior(theta,prior$prior_params))
    }
  }
  
  cores=min(gen_reps,parallel::detectCores())
  if(make_cluster & (gen_parallel | mets_parallel)){
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  mcmc = Metro_Hastings_Stochastic(li_func = llh_mh_adaptive, pars = prior$theta_est, fuel = fuel, prior = prior, prop_sigma = prop.sigma,
                                   par_names = prior$names,
                                   iterations = n.samples, burn_in = n.burn, adapt_par = adapt.par, quiet = !verbose,
                                   ABC = ABC, ABC_eps = ABC_eps, ABC_wt = ABC_wt)
  if(!is.null(fixed)){
    for(i in 1:length(fixed)){
      fixed_idx <- match(names(fixed), prior$names)
      mcmc$trace[,fixed_idx] = fixed[i]
    }
  }
  mcmc$prior = prior
  class(mcmc) = c('list','fuelsgen_mcmc')
  return(mcmc)
}

# ---- rejection ABC driver ----
# generates a bunch of samples from the prior, and evaluates their ABC distance. This is generally a very inefficient way to sample unless the prior looks alot like the posterior--most sampled parameter sets from the prior will have large ABC distance.
abc_rejection <- function(y_obs, fuel, prior,
                          thetas = NULL,
                          N = 1000,
                          eps = NULL,                 # if NULL, set from quantile of distances
                          q_eps = 0.1,               # used only if eps is NULL
                          ABC_wt = c(1, .1),
                          gen_reps = 25,
                          GP.init.size = 32,
                          I.transform = 'exp',
                          mets_parallel = TRUE,
                          make_cluster = TRUE,
                          use_crn = TRUE, crn_base = 4242L,
                          verbose = TRUE)
{
  stopifnot(is.matrix(y_obs$mets))
  
  if(make_cluster & mets_parallel)
    parallel::makeCluster(parallel::detectCores())
  
  # Common-random-numbers for determinism per theta (optional)
  set_seed_from_theta <- function(theta, base = crn_base){
    if (!use_crn) return(invisible(NULL))
    s <- as.integer(abs(round(1e6 * sum(theta))) %% .Machine$integer.max)
    set.seed(base + s)
  }
  
  # 1) draw N prior samples (joint)
  if (verbose) message("ABC-rejection: sampling ", N, " prior draws...")
  if(is.null(thetas))
    thetas <- sample_from_prior(prior,N)
  colnames(thetas) <- prior$names
  
  # --- 2) simulate & distance with 10% status updates ---
  if (verbose) message("Simulating and computing distances (", N, " evaluations)...")
  D <- rep(NA_real_, N)
  marks <- unique(pmax(1, floor(seq(0.1, 1.0, by = 0.1) * N)))
  next_mark_idx <- 1L
  t0 <- Sys.time()
  
  if (verbose) message("Simulating ", N, " prior draws...")
  for (i in seq_len(N)) {
    theta <- thetas[i, ]
    set_seed_from_theta(theta)
    
    sim_fuels <- gen_data(theta,
                          fuel$dimX, fuel$dimY, 1,
                          fuel$X.locs, fuel$X.vals, fuel$Beta,
                          gen_reps, GP.init.size, NULL,
                          I.transform, .217622, FALSE)
    
    mets <- get_mets(sim_fuels, y_obs$info, mets_parallel, make_cluster = FALSE)
    mets <- fix_nan_metrics(theta, mets, fuel)
    
    D[i] <- abc_distance(y_obs$mets, mets$mets, w_met = ABC_wt[1], w_cnt = ABC_wt[2])
    
    # progress ping at each 10%
    if (verbose && next_mark_idx <= length(marks) && i == marks[next_mark_idx]) {
      pct <- 100 * i / N
      dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rate <- i / max(dt, 1e-9)  # iters per sec
      eta  <- (N - i) / max(rate, 1e-9)
      message(sprintf("  ... %3.0f%% (%d/%d), elapsed %.1fs, ETA %.1fs",
                      pct, i, N, dt, eta))
      next_mark_idx <- next_mark_idx + 1L
    }
  }
  
  # 3) choose epsilon if not supplied
  if (is.null(eps)) {
    eps <- as.numeric(quantile(D[is.finite(D)], probs = q_eps))
    if (verbose) message(sprintf("Set eps to %.4g (q_eps=%.2f)", eps, q_eps))
  } else if (verbose) {
    message(sprintf("Using provided eps = %.4g", eps))
  }
  
  # 4) keep acceptances
  keep <- which(is.finite(D) & (D < eps))
  if (verbose) message("Accepted ", length(keep), " / ", N)
  
  list(
    theta_all   = thetas,
    dist_all    = D,
    eps         = eps,
    accepted_ix = keep,
    theta_acc   = if (length(keep)) thetas[keep, , drop = FALSE] else matrix(numeric(0), nrow=0, ncol=ncol(thetas))
  )
}

abc_distance <- function(y_obs_mets, sim_mets,
                         w_met = 1, w_cnt = 1) {
  cmp = c('Exp(n)','Var(n)')
  metric_cols = which(!(colnames(y_obs_mets) %in% cmp))
  cmp_cols = which((colnames(y_obs_mets) %in% cmp))
  
  stopifnot(is.matrix(y_obs_mets), is.matrix(sim_mets))
  Y <- y_obs_mets[, metric_cols, drop=FALSE]
  M <- sim_mets[,   metric_cols, drop=FALSE]
  
  # standardize metrics by observed mean/sd
  muY <- colMeans(Y); sdY <- apply(Y, 2, sd)
  sdY[sdY <= 0 | !is.finite(sdY)] <- 1  # guard
  Yz <- scale(Y, center = muY, scale = sdY)
  Mz <- scale(M, center = muY, scale = sdY)
  
  # energy distance between two samples
  # E = 2 E||X - Y|| - E||X - X'|| - E||Y - Y'||
  # unbiased finite-sample estimator:
  pair_mean <- function(A, B) mean(as.matrix(dist(rbind(A,B)))[1:nrow(A), (nrow(A)+1):(nrow(A)+nrow(B))])
  within_mean <- function(A)  mean(as.matrix(dist(A))[upper.tri(matrix(0,nrow(A),nrow(A)))])
  # efficient versions exist; for moderate sizes this is fine:
  E_xy  <- pair_mean(Yz, Mz)
  E_xx  <- within_mean(Yz)
  E_yy  <- within_mean(Mz)
  E_dist <- 2*E_xy - E_xx - E_yy
  if(!is.finite(E_dist)) E_dist <- .Machine$double.xmax/10
  
  D_total <- w_met * E_dist
  
  if(length(cmp_cols)>0){
    mu_obs = mean(y_obs_mets[,cmp_cols[1]])
    s2_obs = mean(y_obs_mets[,cmp_cols[2]])
    mu_sim = mean(sim_mets[,cmp_cols[1]])
    s2_sim = mean(sim_mets[,cmp_cols[2]])
    
    J <- nrow(Mz)
    # sampling SDs
    sd_mu  <- sqrt(max(s2_sim, 1e-12) / max(J,1))
    sd_s2  <- sqrt(2 / max(J-1,1)) * max(s2_sim, 1e-12)  # chi-square approx
    
    z_mu <- (mu_obs - mu_sim) / ifelse(sd_mu > 0, sd_mu, 1)
    z_s2 <- (s2_obs - s2_sim) / ifelse(sd_s2 > 0, sd_s2, 1)
    
    D_total = D_total + w_cnt * (z_mu^2 + z_s2^2)
  }
  
  return(D_total)
}

fix_nan_metrics = function(theta,metrics,fuel){
  # K function can return NaN's and there is no obvious solution so if it happens we need to re-simulate those processes
  tmp = which(is.na(metrics$mets) | is.nan(metrics$mets) | is.infinite(metrics$mets),arr.ind = T)
  if(length(tmp)>0){
    for(j in 1:nrow(tmp)){
      row = tmp[j,1]
      m = rep(NaN,ncol(metrics$mets))
      while(any(is.nan(m))|any(is.infinite(m))){
        f = gen_data(theta,fuel$dimX,fuel$dimY,1,fuel$X.locs,fuel$X.vals,fuel$Beta,1,fuel$GP.init.size,NULL,fuel$I.transform,.217622,F)
        m = get_mets(f, metrics$info, mets_parallel, make_cluster = F)$mets
      }
      metrics$mets[row,] = m
      # for these, En & Vn are set to n by default because only 1 data point exists, correct the metrics posthoc
      if(info$En){
        metrics$mets[,which(metrics$info$names=='Exp(n)')] = mean(metrics$mets[,which(metrics$info$names=='emp lambda')]*fuel$dimX*fuel$dimY)
      }
      if(info$Vn){
        metrics$mets[,which(metrics$info$names=='Var(n)')] = var(metrics$mets[,which(metrics$info$names=='emp lambda')]*fuel$dimX*fuel$dimY)
      }
    }
  }
  return(metrics)
}
# define prior functions that all take the same inputs
my_dtruncnorm = function(x,params){
  truncdist::dtrunc(x,'norm',a=params[3],b=params[4],mean=params[1],sd=params[2],log=T)
}
my_dunif = function(x,params){
  dunif(x,params[1],params[2],log=T)
}
my_dhalfcauchy = function(x,params){
  LaplacesDemon::dhalfcauchy(x,params,log=T)
}
my_dgamma = function(x,params){
  dgamma(x,params[1],params[2],log=T)
}

lprior = function(theta,prior_params){
  stopifnot("number of parameters not equal to length of parameter info"=length(theta)==length(prior_params))
  lprior = 0
  for(i in 1:length(prior_params)){
    pfunc = switch(prior_params[[i]]$dist,
                   'uniform' = my_dunif,
                   'truncnorm' = my_dtruncnorm,
                   'hcauchy' = my_dhalfcauchy,
                   'gamma' = my_dgamma)
    lprior = lprior + pfunc(theta[i],prior_params[[i]]$params)
  }
  return(lprior)
}

# A minimum of 3 points to define a cluster and points must be within 2m
find_clusters <- function(df, eps=2, minPts=3) {
  # Run DBSCAN
  dbscan_result <- suppressWarnings(dbscan::dbscan(df, eps = eps, MinPts = minPts))
  # get number of clusters
  nclust = length(unique(dbscan_result$cluster)) - 1
  return(nclust)
}

find_connected_components <- function(df) {
  # Function to check if two discs intersect
  discs_intersect <- function(disc1, disc2) {
    # Check if the distance between their centers is less than the sum of their radii
    distance_between_centers <- sqrt((disc1$Xnat - disc2$Xnat)^2 + (disc1$Ynat - disc2$Ynat)^2)
    return (distance_between_centers <= disc1$r + disc2$r)
  }
  
  num_discs <- nrow(df)
  graph <- vector("list", num_discs)
  visited <- rep(FALSE, num_discs)
  components <- 0
  
  # Build the graph
  for (i in 1:(num_discs-1)) {
    for (j in (i+1):num_discs) {
      if(j==i)
        next
      if (discs_intersect(df[i, ], df[j, ])) {
        graph[[i]] <- c(graph[[i]], j)
        graph[[j]] <- c(graph[[j]], i)
      }
    }
  }
  
  # DFS to find connected components
  for (node in 1:num_discs) {
    if (!visited[node]) {
      if(!is.null(graph[[node]])){
        components = components + 1
        visited[node] = T
        visited[graph[[node]]] = T
      }
    }
  }
  
  return (components)
}

moran_geary = function(afs, rq="rook", moran=T, geary=T)
{
  nm = length(afs)
  # previously
  # dd = dim(afs)[1]
  # ids = expand.grid(1:dd, 1:dd)
  # i think this was an error which caused matrix size bugs for rectanglular domains
  dd = dim(afs)
  ids = expand.grid(1:dd[1], 1:dd[2])
  dm = as.matrix(distances::distances(ids))

  if (rq == "rook"){
    wmat = ifelse(dm > 1, 0, 1)
  } else if(rq == "queen"){
    wmat = ifelse(dm > sqrt(2), 0, 1)
  }

  ret = c()
  # Moran's I
  if(moran){
    xminusxbar = c(afs - mean(afs))
    xijm = tidyr::expand_grid(xminusxbar, xminusxbar)
    xtx = matrix((xijm[,1] * xijm[,2])[,1], ncol=nm)
    mi = ( nm * sum(wmat * xtx) ) /
      ( sum(xminusxbar^2) * sum(wmat) )
    ret = c(ret,mi)
  }
  # Geary's C
  if(geary){
    af_vec = c(afs)
    xijg = tidyr::expand_grid(af_vec, af_vec)
    ximinusxj2mat = matrix(((xijg[,1] - xijg[,2])[,1])^2, ncol = nm)
    gc = ((nm - 1) * sum(wmat * ximinusxj2mat)) /
      (2 * sum(wmat) * sum((af_vec - mean(afs))^2))
    ret = c(ret,gc)
  }

  return(ret)
}

# function to remove all shrubs that are completely contained within another shrub
# not removing these rows makes logic for computing segments and breaks significantly more difficult
# this function should only be called from the transect function
remove_rows = function(pp_j,dim)
{
  remove = numeric(nrow(pp_j))
  for(i in 1:nrow(pp_j)){
    for(j in 1:nrow(pp_j)){
      if(i==j)
        next
      if(dim == 'X'){
        if(pp_j$Ypr[i] <= pp_j$Ypr[j] & pp_j$Ymr[i] >= pp_j$Ymr[j]){
          remove[i] = 1
          # we know we are removing row i so break out of j and go to next i
          break
        }
      } else{
        if(pp_j$Xpr[i] <= pp_j$Xpr[j] & pp_j$Xmr[i] >= pp_j$Xmr[j]){
          remove[i] = 1
          # we know we are removing row i so break out of j and go to next i
          break
        }
      }
    }
  }
  return(pp_j[!remove,])
}
transect = function(t_locs,dims,dat,dimX,dimY)
{
  n_transects = length(t_locs)
  segment_lengths = vector(mode='list',n_transects)
  break_lengths = vector(mode='list',n_transects)

  for(j in 1:n_transects){
    if(dims[j]=='X'){
      x_dist_from_transect = abs(dat$Xnat - t_locs[j])
      x_within_r_from_transect = x_dist_from_transect <= dat$r

      # calculate segment lengths
      pp_j = dat[x_within_r_from_transect,]
      if(nrow(pp_j)>0){
        # pmax(0,pp_j$Y - pp_j$r) takes care of the case where the left boundary is < 0
        # pmin(pp_j$Y + pp_j$r,dimY) takes care of the case where the right boundary is > dimY
        pp_j$Ymr = pmax(0,pp_j$Y - pp_j$r)
        pp_j$Ypr = pmin(pp_j$Y + pp_j$r,dimY)
        # sort shrubs in order of left boundary
        pp_j = pp_j[order(pp_j$Ymr, pp_j$Ypr),]
        # remove shrubs that are completely contained within previous shrub
        pp_j = remove_rows(pp_j,dims[j])
        segment_lengths[[j]] = pp_j$Ypr[1] - pp_j$Ymr[1]
        if(pp_j$Ymr[1]>0){
          break_lengths[[j]] = c(pp_j$Ymr[1])
        } else{
          break_lengths[[j]] = numeric()
        }
        if(nrow(pp_j)>1){
          k = 1
          for(i in 2:nrow(pp_j)){
            # this is an oversimplification that does not account for segment length of the line bisecting the circle
            if(pp_j$Ypr[i-1] == dimY){
              # end if we reached the boundary
              break
            }
            # check if circle is overlapping previous
            if( pp_j$Ymr[i] <= pp_j$Ypr[i-1]){
              # update segment lengths to reflect total length of these connected circles
              segment_lengths[[j]][k] = segment_lengths[[j]][k] + (pp_j$Ypr[i] - pp_j$Ypr[i-1])
            } else if(pp_j$Ymr[i] > pp_j$Ypr[i-1]){
              segment_lengths[[j]] = c(segment_lengths[[j]], pp_j$Ypr[i] - pp_j$Ymr[i])
              break_lengths[[j]] = c(break_lengths[[j]],pp_j$Ymr[i] - pp_j$Ypr[i-1])
              k = k + 1
            }
          }
        }
        d_to_bound = dimY - pp_j$Ypr[nrow(pp_j)]
        if(d_to_bound>0){
          break_lengths[[j]] = c(break_lengths[[j]],d_to_bound)
        }
      } else{
        # no shrubs touch transect a.k.a. one break of length dimY and no segments
        segment_lengths[[j]] = 0
        break_lengths[[j]] = dimY
      }

    } else if(dims[j]=='Y'){
      y_dist_from_transect = abs(dat$Ynat - t_locs[j])
      y_within_r_from_transect = y_dist_from_transect <= dat$r

      # calculate segment lengths
      pp_j = dat[y_within_r_from_transect,]
      if(nrow(pp_j)>0){
        pp_j$Xmr = pmax(0,pp_j$X - pp_j$r)
        pp_j$Xpr = pmin(pp_j$X + pp_j$r,dimX)
        pp_j = pp_j[order(pp_j$Xmr),]
        # remove shrubs that are completely contained within previous shrub
        pp_j = remove_rows(pp_j,dims[j])
        segment_lengths[[j]] = pp_j$Xpr[1] - pp_j$Xmr[1]
        if(pp_j$Xmr[1]>0){
          break_lengths[[j]] = c(pp_j$Xmr[1])
        } else{
          break_lengths[[j]] = numeric()
        }
        if(nrow(pp_j)>1){
          k = 1
          for(i in 2:nrow(pp_j)){
            # this is an oversimplification that does not account for segment length of the line bisecting the circle
            if(pp_j$Xpr[i-1] == dimX){
              break
            }
            if(pp_j$Xpr[i] <= pp_j$Xpr[i-1]){
              # move to next as this shrub is completely contained in the previous shrub
              next
            }
            # check if cirle is overlapping previous
            if( pp_j$Xmr[i] <= pp_j$Xpr[i-1] ){
              # update segment lengths to reflect total length of these connected circles
              segment_lengths[[j]][k] = segment_lengths[[j]][k] + (pp_j$Xpr[i] - pp_j$Xpr[i-1])
            } else{
              segment_lengths[[j]] = c(segment_lengths[[j]],pp_j$Xpr[i] - pp_j$Xmr[i])
              break_lengths[[j]] = c(break_lengths[[j]],pp_j$Xmr[i] - pp_j$Xpr[i-1])
              k = k + 1
            }
          }
        }
        d_to_bound = dimX - pp_j$Xpr[nrow(pp_j)]
        if(d_to_bound>0){
          break_lengths[[j]] = c(break_lengths[[j]],d_to_bound)
        }
      } else{
        segment_lengths[[j]] = 0
        break_lengths[[j]] = dimX
      }
    } else{
      stop("dim must be one of 'X' or 'Y'")
    }
  }
  for(i in 1:n_transects){
    seg_sum = ifelse(length(segment_lengths[[i]])>0, sum(segment_lengths[[i]]), 0)
    break_sum = ifelse(length(break_lengths[[i]])>0, sum(break_lengths[[i]]), 0)
    max = ifelse(dims[i]=='X',dimY,dimX)
    if(abs(seg_sum + break_sum - max)>1e-4){
      #cat('Location: ',t_locs[i],'\nDim: ',dims[i],'\nSegments: ',round(segment_lengths[[i]],2),'\nBreaks: ',round(break_lengths[[i]],2),'\nSum: ',seg_sum+break_sum,'\n')
      warning('Error in transect ',i,': sum of segments and breaks != ',max,'\n')
    }
  }
  return(list(segments = segment_lengths,
              breaks = break_lengths))
}

#' @title Compute metrics for fuels
#'
#' @description Compute metrics for each fuel realization in fuels object
#' @param fuels fuel object outputs from gen_fuels()
#' @param info metrics info from get_mets_info(), if NULL, default metrics are used
#' @param parallel compute metrics in parallel
#' @param make_cluster make parallel cluster in function
#' @export
#'
get_mets = function(fuels, info=NULL, parallel = F, make_cluster = F)
{
  if(is.null(info)){
    info = get_mets_info()
  }
  if(fuels$reps>1 & parallel){
    metrics = mets_parallel(fuels, info, make_cluster)
  } else{
    metrics = c()
    for(i in 1:fuels$reps){
      metrics = rbind(metrics,mets(fuels$dat[[i]], fuels$W, info))
    }
  }
  if(info$En){
    En = jitter(rep(ifelse(fuels$reps>1,mean(metrics[,which(info$names=='emp lambda')]*fuels$dimX*fuels$dimY,na.rm=T),metrics[,which(info$names=='emp lambda')]*fuels$dimX*fuels$dimY),nrow(metrics)),amount = .001)
  }
  if(info$Vn){
    Vn = jitter(rep(ifelse(fuels$reps>1,var(metrics[,which(info$names=='emp lambda')]*fuels$dimX*fuels$dimY,na.rm=T),metrics[,which(info$names=='emp lambda')]*fuels$dimX*fuels$dimY),nrow(metrics)),amount = .001)
  }
  if(info$En){
    metrics = cbind(metrics, En)
  }
  if(info$Vn){
    metrics = cbind(metrics, Vn)
  }
  colnames(metrics) = c(info$names,rep("",ncol(metrics)-length(info$names)))
  
  return(list(mets=metrics,info=info))
}

mets_parallel = function(fuels,info,make_cluster)
{
  if(make_cluster){
    cl = parallel::makeCluster(min(fuels$reps,parallel::detectCores()))
    doParallel::registerDoParallel(cl)
  }
  metrics = foreach::foreach(i=1:fuels$reps,.combine = 'rbind') %dopar% mets(fuels$dat[[i]],fuels$W, info)
  if(make_cluster)
    parallel::stopCluster(cl)
  return(metrics)
}

mets = function(dat, W, info)
{
  dimX = W$dimX
  dimY = W$dimY
  matsplitter=function(M, r, c)
  {
    rg = (row(M)-1)%/%r+1
    cg = (col(M)-1)%/%c+1
    rci = (rg-1)*max(cg) + cg
    N = prod(dim(M))/r/c
    cv = unlist(lapply(1:N, function(x) M[rci==x]))
    dim(cv)=c(r,c,N)
    cv
  }
  
  metrics = c()
  n = nrow(dat)
  if (n == 0){
    remove = which(info$names %in% c("Exp(n)","Var(n)"))
    return(rep(0,length(info$names[-remove])))
  }
  
  #--- NCC ---#
  # I don't think this is working like we want it to. It should be separating better on rho
  if(info$ncc){
    if(n==1){
      ncc.val = 0
      ncc.time = 0
    } else{
      ptm = proc.time()[3]
      # r_mat =  matrix(rep(dat$r, n), ncol=n, byrow=T)
      # r_pairs = r_mat + t(r_mat)
      # cent_dist = as.matrix(distances::distances(dat[, c("Xnat", "Ynat")]))
      # incidence = cent_dist <= r_pairs
      # ig = igraph::graph_from_adjacency_matrix(incidence * (1 / cent_dist),
      #                                          weighted = T, diag = F, mode="undirected")
      #ds = igraph::distances(ig)
      #ds[is.infinite(ds)] = 0
      #max_d = max(ds)
      # ncc.val = igraph::components(ig)$no
      ncc.val = find_connected_components(dat)
      ncc.time = proc.time()[3] - ptm
    }
  }
  #--- N Clusters ---#
  if(info$nclust){
    if(n==1){
      nclust.val = 0
      nclust.time = 0
    } else{
      ptm = proc.time()[3]
      nclust.val = find_clusters(dat[, c("Xnat", "Ynat")])
      nclust.time = proc.time()[3] - ptm
    }
  }
  #--- N holes ---#
  if(info$nholes){
    if(n==1){
      nholes.val = 0
      nholes.time = 0
    } else{
      ptm = proc.time()[3]
      # if(!info$ncc){
      # these used to be computed in ncc so if we don't do ncc we must compute them here
      r_mat =  matrix(rep(dat$r, n), ncol=n, byrow=T)
      r_pairs = r_mat + t(r_mat)
      cent_dist = as.matrix(distances::distances(dat[, c("X", "Y")]))
      # }
      
      # why a threshold of 5?
      homo = TDAstats::calculate_homology(cent_dist / r_pairs,
                                          format = "distmat", threshold = 5, return_df = T)
      #homo = TDApplied::PyH(cent_dist / r_pairs,
      #                      distance_mat = TRUE, thresh = 5, ripser=ripser)
      # nholes.val = nrow(dplyr::filter(homo, dimension==1 & birth < 1 & death > 1))
      nholes.val = nrow(dplyr::filter(homo, dimension==1))
      nholes.time = proc.time()[3] - ptm
    }
  }
  #--- Gridded area ---#
  if(info$grid.area | info$moran | info$geary){
    ptm = proc.time()[3]
    dx = dimX/10
    dy = dimY/10
    stepsize = dx/4
    xstart = 0

    area_fractions = matrix(ncol=dimX / dx, nrow=dimY / dy)

    for (i in 1:(dimX / dx)) {
      ystart = 0
      for (j in 1:(dimY / dx)) {
        my_grid = pracma::meshgrid(seq(xstart, xstart + dx - stepsize, by=stepsize),
                                   seq(ystart, ystart + dy - stepsize, by=stepsize))
        mc_d = sqrt(plgp::distance(cbind(c(my_grid$X), c(my_grid$Y)), dat[,c("Xnat", "Ynat")]))

        area_fractions[j,i] = mean(apply(mc_d < dat$r, 1, any))
        ystart = ystart + dy
      }
      xstart = xstart + dx
    }
    af.time = proc.time()[3] - ptm
    total_area = mean(area_fractions) * dimX * dimY
    empty_cells = sum(area_fractions == 0)
    full_cells = sum(area_fractions == 1)
    af_var = var(c(area_fractions))
  }

  #--- Moran I ---#
  if(info$moran | info$geary){
    ptm = proc.time()[3]
    af0 = ifelse(area_fractions == 0, 0, 1)
    if(info$rook){
      r1 = moran_geary(area_fractions, "rook", info$moran, info$geary)
      r01 = moran_geary(af0, "rook", info$moran, info$geary)
    }
    if(info$queen){
      q1 = moran_geary(area_fractions, "queen", info$moran, info$geary)
      q01 = moran_geary(af0, "queen", info$moran, info$geary)
    }

    af_half = matrix(apply(matsplitter(area_fractions, 2, 2), 3, mean),
                     ncol = dimX / (dx * 2), byrow = T)
    af_half0 = ifelse(af_half == 0, 0, 1)

    if(info$rook){
      r2 = moran_geary(af_half, "rook", info$moran, info$geary)
      r02 = moran_geary(af_half0, "rook", info$moran, info$geary)
    }
    if(info$queen){
      q2 = moran_geary(af_half, "queen", info$moran, info$geary)
      q02 = moran_geary(af_half0, "queen", info$moran, info$geary)
    }

    # what does this mean? These are 2-vectors so x[is.nan(x)] could be length 1 or length 2
    # If it's length 1, we are setting it to zero always and if its length 2 we set to c(0,1)
    # so this is saying in words "if 1 of them is NaN, make that one zero, and if both are Nan, set the moran to 0 and the geary to 1"
    # Need to ask Braden since it's not intuative. Try with commented out to see if NaN's are a problem,.
    #
    # this  might be saying set moran NaN's to 0 and Geary NaN's to 1
    if(info$rook){
      r01[is.nan(r01)] = c(0, 1)
      r02[is.nan(r02)] = c(0, 1)
    }
    if(info$queen){
      q1[is.nan(q1)] = 1
      q2[is.nan(q2)] = 1
      q01[is.nan(q01)] = 1
      q02[is.nan(q02)] = 1
    }
    moran.time = proc.time()[3] - ptm
  }
  #--- Perimeter ---#
  if(info$perim){
    ptm = proc.time()[3]
    r_sum = sum(dat$r)
    store = c(0, 0)

    for (i in 1:n){
      n_perm = max(c(1, round((dat$r[i] / r_sum) * 5000)))
      rads = seq(0, 2 * pi, length.out=(n_perm + 1))[1:n_perm]
      xy = rep(c(dat[i,1], dat[i,2]), each=n_perm)
      this = cbind(dat$r[i] * cos(rads) + dat[i,1], dat$r[i] * sin(rads) + dat[i,2])
      dis = sqrt(plgp::distance(this, dat[,c("Xnat", "Ynat")]))
      perim_pts = sum(abs(matrixStats::rowMins(dis) - dat$r[i]) < 0.000001)

      store[1] = store[1] + n_perm
      store[2] = store[2] + perim_pts
    }

    pct_perm = store[2] / store[1]
    perimeter = pct_perm * sum(2 * dat$r * pi)
    perim.time = proc.time()[3] - ptm
  }
  #--- Transects ---#
  if(info$transect){
    ptm = proc.time()[3]
    if(is.null(info$dim.trans)){
      if(info$n.trans==1){
        info$dim.trans = 'X'
      } else if(info$n.trans>1){
        if(info$n.trans%%2==0){
          info$dim.trans = c(rep('X',as.integer(info$n.trans/2)),rep('Y',as.integer(info$n.trans/2)))
        } else{
          info$dim.trans = c(rep('X',as.integer(info$n.trans/2)+1),rep('Y',as.integer(info$n.trans/2)))
        }
      } else{
        stop('n.trans<1')
      }
    }
    if(is.null(info$loc.trans)){
      # evenly space the transects in each dimension
      if(info$n.trans==1){
        info$loc.trans = dimX/2
      } else{
        if(info$n.trans%%2==0){
          n.x = info$n.trans/2; n.y = info$n.trans/2
        } else{
          n.x = as.integer(info$n.trans/2) + 1; n.y = as.integer(info$n.trans/2)
        }
        x_start = dimX/(info$n.trans/2+1); x_end = dimX - x_start
        y_start = dimY/(info$n.trans/2+1); y_end = dimX - x_start
        info$loc.trans = c(seq(x_start,x_end,length.out=n.x),seq(y_start,y_end,length.out=n.y))
      }
    }
    t = transect(info$loc.trans,info$dim.trans,dat,dimX,dimY)
    segments = unlist(t$segments)
    breaks = unlist(t$breaks)
    transect.time = proc.time()[3] - ptm
  }

  if(info$ncc){
    metrics = c(metrics,ncc.val)
  }
  if(info$nclust){
    metrics = c(metrics,nclust.val)
  }
  if(info$nholes){
    metrics = c(metrics,nholes.val)
  }
  if(info$grid.area){
    metrics = c(metrics, total_area, full_cells, empty_cells, af_var * 100)
  }
  if(info$moran | info$geary){
    if(info$rook){
      metrics = c(metrics,c(r1,r2,r01,r02))
    }
    if(info$queen){
      metrics = c(metrics,c(q1,q2,q01,q02))
    }
  }
  if(info$perim){
    metrics = c(metrics,perimeter)
  }
  if(info$transect){
    nseg = length(segments)
    nbreak = length(breaks)
    metrics = c(metrics,nseg,ifelse(nseg>0,mean(segments),0),ifelse(nseg>1,var(segments),0),
                     nbreak,ifelse(nbreak>0,mean(breaks),0),ifelse(nbreak>1,var(breaks),0))
  }

  # Summary statistics 
  if(info$Kinhom | info$pcf | info$Ginhom){
    # something here is calling glm() which is ocationally throwing an error
    P = spatstat.geom::ppp(x = dat$X,y = dat$Y, window = W)
    ppm = spatstat.model::ppm # this is a strange workaround i had to do for the mcmc function to find the ppm function, no idea why
    if(P$n>10){
      trend = ppm(P ~ polynom(x,y,3), interaction = spatstat.model::Poisson())
      lambda = predict(trend, type='trend')
    } else if(P$n>5){
      trend = ppm(P ~ polynom(x,y,2), interaction = spatstat.model::Poisson())
      lambda = predict(trend, type='trend')
    } else if(P$n>2){
      trend = ppm(P ~ polynom(x,y,1), interaction = spatstat.model::Poisson())
      lambda = predict(trend, type='trend')
    } else{
      lambda = NULL
    }
    if(info$Kinhom | info$pcf){
      Kint = tryCatch({
        Kin = spatstat.explore::Kinhom(P,window = W,lambda = NULL,
                                       correction = 'border', normpower = 2)
        Kin$border = spatstat.model::safePositiveValue(Kin$border)
        # is r_max always the same for the same window? Check this, otherwise the integrals won't always match
        suppressWarnings({spatstat.univar::integral(Kin)[2]^(1/4)})
      }, error = function(e) {
        0
      })
    }
    if(info$Kinhom){
      metrics = c(metrics,Kint)
    }
    if(info$pcf){
      Pcf = spatstat.explore::pcf.fv(Kin)
      Pcf$pcf = spatstat.model::safePositiveValue(Pcf$pcf)
      PCFint = suppressWarnings({spatstat.univar::integral(Pcf^(1/4))[2]})
      metrics = c(metrics,PCFint)
    }
    if(info$Ginhom){
      G = spatstat.explore::Ginhom(P,lambda = lambda)
      G$bord = spatstat.model::safePositiveValue(G$bord)
      Gint = suppressWarnings({spatstat.univar::integral(G^(1/4))[2]})
      metrics = c(metrics,Gint)
    }
  }
  
  # include empirical estimates of parameters as metrics
  metrics = c(metrics,n/dimX/dimY,mean(dat$r),ifelse(n==1,0,var(dat$r)))
  metrics = jitter(metrics,amount = .001) # zeros cause issue with cov matrix
  return(metrics)
}

#' @title Setup prior for MCMC
#'
#' @description Compute necessary prior information for MCMC sampling
#' @param fuel observed fuel object, output from gen_fuels()
#' @param metrics metrics for fuel object, output from get_mets()
#' @param est_cov_obs estimate covariance using only observations (not recommended for n<10)
#' @param est_cov_samples number of samples from priors used for augmenting observed metrics, only used if est_cov_obs=F
#' @param est_cov_reps number of fuel replicates to generate at each sample, only used if est_cov_obs=F
#' @param est_rho_prior list containing string 'prior' denoting the prior distribution and vector 'params' containing the parameters of the distribution for the lengthscale. Currently only 'unif' and 'gamma' priors are supported.
#' @param gen_parallel generate the fuels in parallel (can be faster for large domains with many fuel elements, but will be slower for small domains)
#' @param mets_parallel compute metrics in parallel (almost always faster)
#' @param make_cluster should a parallel cluster be created within the function?
#' @param seed random seed for reproducible fuel generation and metrics
#' @returns List containing necessary prior precomputing for MCMC
#' @export
#'
get_prior_info = function(fuel,metrics,est_cov_obs=T,sigma_diag=F,
                          est_cov_samples=25,est_cov_reps=25,
                          est_rho_prior=list(prior='gamma',params=c(1,10/fuel$dimX)),
                          mu_mean = 1.5, mu_sd = .5,
                          gen_parallel = F, mets_parallel = T, make_cluster = F,
                          eps = F,
                          I.transform='logistic',GP.init.size = 64, seed = NULL)
{
  cluster.made = 0
  if(!est_cov_obs & (gen_parallel | mets_parallel)){
    cores = min(max(est_cov_samples,est_cov_reps),parallel::detectCores())
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    cluster.made = 1
  }
  lambda_est = numeric(fuel$reps)
  mu_est = numeric(fuel$reps)
  sigma_est = numeric(fuel$reps)
  for(i in 1:fuel$reps){
    lambda_est[i] = nrow(fuel$dat[[i]])/fuel$dimX/fuel$dimY
    mu_est[i] = mean(fuel$dat[[i]]$r)
    sigma_est[i] = sd(fuel$dat[[i]]$r)
  }
  lambda_min = min(lambda_est)
  lambda_max = max(lambda_est)
  dispersion_est = mean(lambda_est*fuel$dimX*fuel$dimY)/var(lambda_est*fuel$dimX*fuel$dimY) # nu
  lambda_est = mean(lambda_est)
  
  if(fuel$reps>=5){
    # mean paramaterized cmp model
    counts = numeric()
    for(i in 1:fuel$reps){counts[i] = nrow(fuel$dat[[i]])}
    mpcmp.mod = mpcmp::glm.cmp(counts~1)
    
    com_mu_est = mpcmp.mod$lambda[1]
    dispersion_est = mpcmp.mod$nu
  } else{
    # fix dispersion at 1
    dispersion_est = 1
  }
  
  mu_est = mean(mu_est,na.rm=T) # there may be observed domains with 0 fuel elements
  sigma_est = mean(sigma_est,na.rm=T) # there may be observed domains with 0 fuel elements
  # rho_max = max(fuel$dimX,fuel$dimY)/2
  rho_max = 1 # learning on standard domain
  
  # gamma parameters s.t. q95 ~~ rho_max
  #rho_prior = "gamma"
  rho_prior = "uniform"
  rho_shape = NULL # 1
  rho_rate = NULL # .3

  # set bounds for each parameters
  
  # the lower bound on lamda is very important. Certain metrics calculations will fail if not enough trees are sampled.
  # The error we get with small lamda is "Point cloud must have at least 2 points and at least 2 dimensions." I need to find out
  # which metric calculation is causing this.
  # We choose to assume that the true lambda will not be smaller than 1/2 the observed relative density or greater than 2 times the max observed relative density
  
  if(est_cov_obs){
    if(fuel$reps>1){
      # estimate Sigma matrix using observed metrics
      Sigma = cov(metrics$mets)
      nu = nrow(metrics$mets)
      if(sigma_diag)
        Sigma = diag(diag(Sigma))
    } else{
      # prior covariance is 10% standard error
      Sigma = diag(diag(as.numeric(.1^2*metrics$mets)))
      nu = 2 # this isn't quite true, but we need nu>1 for the prior (nu-1)*Sigma to work
    }
    metrics_samples = c()
  } else{
    if(!is.null(seed))
      set.seed(seed)
    # use lambda_est, mu_est, sigma_est + prior samples of rho to estimate Sigma
    lambda_samples = truncdist::rtrunc(est_cov_samples,'norm',a=0,mean=lambda_est,sd=.1*lambda_est)
    mu_samples = truncdist::rtrunc(est_cov_samples,'norm',a=0,b=3,mean=mu_est,sd=.1*mu_est)
    sigma_samples = truncdist::rtrunc(est_cov_samples,'norm',a=0,mean=sigma_est,sd=.1*sigma_est)
    if(est_rho_prior$prior == 'gamma'){
      rho_samples = truncdist::rtrunc(est_cov_samples,'gamma',a=0,b=rho_max,
                                      shape=est_rho_prior$params[1],rate=est_rho_prior$params[2])
    } else if(est_rho_prior$prior == 'unif'){
      rho_samples = runif(est_cov_samples,est_rho_prior$params[1],est_rho_prior$params[2])
    } else{
      stop('enter either gamma or unif as rho prior')
    }
    dispersion_samples = truncdist::rtrunc(est_cov_samples,'norm',a=0,mean=dispersion_est,sd=.1*dispersion_est)
    metrics_samples = c()
    for(i in 1:est_cov_samples){
      # generate new fuel realizations using sampled parameters and same inputs as observations if available
      theta = c(rho_samples[i],mu_samples[i],sigma_samples[i],lambda_samples[i],dispersion_samples[i])
      fuel_samp = gen_data(theta,
                           fuel$dimX, fuel$dimY, fuel$heterogeneity.scale, 
                           fuel$X.locs, fuel$X.vals, fuel$Beta, 
                           est_cov_reps, GP.init.size, seed = seed, I.transform, fuel$logis.scale, gen_parallel)
      M = get_mets(fuel_samp,metrics$info,mets_parallel, make_cluster = F)$mets
      tmp = which(is.nan(M) | is.infinite(M),arr.ind = T)
      if(length(tmp)>0){
        for(j in 1:nrow(tmp)){
          row = tmp[j,1]
          m = rep(NaN,ncol(M))
          while(any(is.nan(m))|any(is.infinite(m))){
            f = gen_data(theta,fuel$dimX,fuel$dimY,1,fuel$X.locs,fuel$X.vals,fuel$Beta,1,GP.init.size,NULL,I.transform,fuel$logis.scale,F)
            m = get_mets(f, metrics$info, mets_parallel, make_cluster = F)$mets
          }
          M[row,] = m
          # for these, En & Vn are set to n by default because only 1 data point exists, correct the metrics posthoc
          if(info$En){
            M[,which(colnames(M)=='Exp(n)')] = mean(M[,which(colnames(M)=='emp lambda')]*fuel$dimX*fuel$dimY)
          }
          if(info$Vn){
            M[,which(colnames(M)=='Var(n)')] = var(M[,which(colnames(M)=='emp lambda')]*fuel$dimX*fuel$dimY)
          }
        }
      }
      # add mean metrics to matrix
      metrics_samples = rbind(metrics_samples,colMeans(M))
    }
    Sigma = cov(rbind(metrics$mets,metrics_samples))
    nu = est_cov_samples * est_cov_reps + nrow(metrics$mets)
    if(diag_sigma)
      Sigma = diag(diag(Sigma))
  }
  
  # set rho_est to prior mean
  if(rho_prior == 'gamma'){
    rho_est = rho_shape/rho_rate
  } else if(rho_prior == 'uniform'){
    rho_est = rho_max/2
  }
  
  if(cluster.made)
    parallel::stopCluster(cl)
  
  Sigma = Sigma + 1e-8*diag(dim(Sigma)[1])
  Sinv = chol2inv(chol(Sigma))
  ldetS = determinant(Sigma)$modulus
  
  # save everything needed to calculate priors
  # prior_params = list()
  # prior_params$lb = prior.lb; prior_params$ub = prior.ub 
  # prior_params$rho_prior = rho_prior
  # prior_params$rho_max = rho_max
  # prior_params$rho_shape = rho_shape
  # prior_params$rho_rate = rho_rate
  # prior_params$mu_prior = "truncnorm"
  # prior_params$mu_mean = mu_mean
  # prior_params$mu_sd = mu_sd
  # prior_params$prec_prior = 'gamma'
  # prior_params$prec_shape = 1
  # prior_params$prec_rate = 1e-3
  # prior_params$var_prior = 'hcauchy'
  # prior_params$var_scale = .1
  # prior_params$lambda_prior = 'gamma'
  # prior_params$lambda_shape = 1
  # prior_params$lambda_rate = 1/lambda_est
  # prior_params$dispersion_shape = dispersion_est
  # prior_params$dispersion_rate = 1 
  # prior_params$dist = c(rho_prior,'truncnorm','hcauchy','gamma','gamma')
  
  if(rho_prior=='uniform'){
    rho_params = c(0,rho_max)
  } else{
    rho_params = c(rho_shape,rho_rate)
  }
  
  names =c('rho',    'mu',       's2',     'lambda','nu')
  dist = c(rho_prior,'truncnorm','hcauchy','gamma', 'gamma')
  theta_est = c(rho_est,mu_est,sigma_est^2,lambda_est,dispersion_est)
  lb = c(0,       0, 0,              lambda_min/2,  0.1)
  ub = c(rho_max, 3, 10*sigma_est^2, lambda_max*2,  8)
  
  if(eps){
    names = c(names,'eps')
    dist = c(dist,'hcauchy')
    theta_est = c(theta_est,0)
    lb = c(lb, 0)
    ub = c(ub, 1)
  }
  names(theta_est) = names
  # all info needed to compute prior
  prior_params = list(
    list(name='rho',   dist=rho_prior,  params=rho_params,          bounds=c(0,rho_max)),
    list(name='mu',    dist='truncnorm',params=c(mu_est,2*sigma_est,0,3),bounds=c(0,3)),
    list(name='s2',    dist='hcauchy',  params=c(.1),               bounds=c(0,10*sigma_est^2)),
    list(name='lambda',dist='gamma',
         params=c(2,(2-1)/lambda_est),bounds=c(lambda_min/2,2*lambda_max)), # mode at estimate
    list(name='nu',    dist='gamma',
         params=c(2,1),bounds=c(.1,8))) # mode at poisson
  if(eps){
    prior_params = append(prior_params,list(name='eps',   dist='hcauchy',  params=c(.1), bounds=c(0,1)))
  }
  ret = list('names' = names,
             'dist' = dist,
             'lb' = lb,
             'ub' = ub,
             'prior_params' = prior_params,
             'theta_est' = theta_est,
             'sim_metrics' = metrics_samples,
             'Sigma' = Sigma,
             'nu' = nu, # numer of 'observations' used to estimate covariance
             'Sinv' = Sinv,
             'ldetS' = ldetS)
  class(ret) = c('list','fuelsgen_prior')
  return(ret)
}

#' @title Metrics info
#'
#' @description Get a list of metrics indicators
#' @param ncc 
#' @param nholes number of holes in the binary map
#' @param grid.area 
#' @param moran Moran's I spatial autocorrelation. Moran=F by default because it generally doesn't add extra information from Geary.
#' @param geary Geary's C spacial autocorrelation
#' @param rook rook method for connected components in Moran/Geary. Rook is turned off my default as it generally doesn't add extra information from queen.
#' @param queen queen method for connected components in Moran/Geary.
#' @param perim perimeter of disks. Perimeter=F by default as it's very expensive to compute.
#' @param transect Transect lines
#' @param n.trans Number of transect lines layed out in a grid
#' @details returns a list of indicators for each metric
#' @export
get_mets_info = function(dimX,dimY,
                         ncc=T, nclust=T, nholes=T,
                         grid.area=T,moran=F,geary=T,rook=F,queen=T,
                         perim=F,transect=F,n.trans=0,
                         Kinhom=F,pcf=F,Ginhom=F,
                         En=T,Vn=T)
{
  names = c()
  if(ncc){
    names = c(names,'n cc')
  }
  if(nclust){
    names = c(names, 'n clust')
  }
  if(nholes){
    names = c(names,'n holes')
  }
  if(grid.area){
    names = c(names,'total area','full cells','empty cells','area frac var')
  }
  if(moran | geary){
    if(rook){
      names = c(names,'Moran Geary r1','Moran Geary r2','Moran Geary r01','Moran Geary r02')
    }
    if(queen){
      names = c(names,'Moran Geary q1','Moran Geary q2','Moran Geary q01','Moran Geary q02')
    }
  }
  if(perim){
    names = c(names,'perimeter')
  }
  if(transect){
    names = c(names,'n segs','mean seg length','var seg length',
              'n breaks','mean break length','var break length')
  }
  if(Kinhom){
    names = c(names,'Kinhom_int')
  }
  if(pcf){
    names = c(names,'PCF_int')
  }
  if(Ginhom){
    names = c(names,'Ginhom_int')
  }
  # include empirical estimates of parameters as metrics
  names = c(names,'emp lambda','emp r mean','emp r var')
  
  if(En){
    names = c(names,'Exp(n)')
  }
  if(Vn){
    names = c(names,'Var(n)')
  }
  
  return(list('ncc'=ncc,
              'nclust'=nclust,
              'nholes'=nholes,
              'grid.area'=grid.area,
              'moran'=moran,
              'geary'=geary,
              'rook'=rook,
              'queen'=queen,
              'perim'=perim,
              'transect'=transect,
              'n.trans'=n.trans,
              'Kinhom'=Kinhom,
              'pcf'=pcf,
              'Ginhom'=Ginhom,
              'En'=En,'Vn'=Vn,
              'names'=names))
}

# ---- sample one vector theta from the prior ----
sample_from_prior <- function(prior,N) {
  # ---- helper: robust truncated sampler via truncdist ----
  .rtrunc_safe <- function(n, spec, a, b, ...) {
    if (!is.finite(a)) a <- -Inf
    if (!is.finite(b)) b <-  Inf
    if (!(a < b)) stop("Bad truncation: a>=b (a=", a, ", b=", b, ")")
    truncdist::rtrunc(n, spec, a = a, b = b, ...)
  }
  
  p <- length(prior$names)
  theta <- matrix(nrow=N,ncol=p)
  
  for (i in seq_len(p)) {
    nm   <- prior$names[i]
    dist <- if (!is.null(prior$dist)) prior$dist[i] else "uniform"
    info <- prior$prior_params[[i]]
    bnds <- info$bounds
    pars <- info$params
    
    # Normalize dist labels a bit
    dist <- tolower(dist)
    if (dist %in% c("unif","uniform")) dist <- "uniform"
    if (dist %in% c("tnorm","truncnorm","truncatednormal","trunc_normal")) dist <- "truncnorm"
    
    theta[,i] <- switch(
      dist,
      
      # Uniform on bounds
      "uniform" = runif(N, min = bnds[1], max = bnds[2]),
      
      # Truncated Normal: params = c(mean, sd, a, b) (if a/b missing, fall back to bounds)
      "truncnorm" = {
        a <- if (length(pars) >= 3 && is.finite(pars[3])) pars[3] else bnds[1]
        b <- if (length(pars) >= 4 && is.finite(pars[4])) pars[4] else bnds[2]
        .rtrunc_safe(N, "norm", a = a, b = b, mean = pars[1], sd = pars[2])
      },
      
      # Gamma: params = c(shape, scale); respect bounds via truncation
      "gamma" = .rtrunc_safe(N, "gamma", a = bnds[1], b = bnds[2],
                             shape = pars[1], scale = pars[2]),
      
      # Lognormal: params = c(meanlog, sdlog); truncated to bounds
      "lognormal" = .rtrunc_safe(N, "lnorm", a = bnds[1], b = bnds[2],
                                 meanlog = pars[1], sdlog = pars[2]),
      
      # Beta: params = c(alpha, beta); truncated to bounds (typically [0,1])
      "beta" = .rtrunc_safe(N, "beta", a = bnds[1], b = bnds[2],
                            shape1 = pars[1], shape2 = pars[2]),
      
      # Normal (untruncated): params = c(mean, sd) but clipped to bounds afterward
      "normal" = {
        x <- rnorm(N, mean = pars[1], sd = pars[2])
        min(max(x, bnds[1]), bnds[2])
      },
      
      # hcauchy: params = scale
      "hcauchy" = {
        u = runif(N,LaplacesDemon::phalfcauchy(bnds[1], scale = pars[1]),LaplacesDemon::phalfcauchy(bnds[2], scale = pars[1]))
        LaplacesDemon::qhalfcauchy(u, scale = pars[1])
      },
      
      # default fallback → uniform on bounds
      {
        warning("Unknown prior dist '", dist, "' for ", nm, "; using Uniform(bounds).")
        runif(1, min = bnds[1], max = bnds[2])
      }
    )
  }
  names(theta) <- prior$names
  theta
}
