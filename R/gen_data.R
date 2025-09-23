relu = function(x,m=0,bound=0){
  x[x<=bound] = m
  return(x)
}
plogis = function(x,location,scale){
  x = x - mean(x)
  1 / (1 + exp(-(x-location)/scale))
}

# post-hoc repulsion operator for generated fuels map
# try weighting the repulsion effect by the radius
# If boundary = 'periodic', periodic boundary conditions are used
# to calculate distance. This helps prevent points from pushing up
# against the boundary, which is bad for tiling.
repulsion = function(fuels, eps, B, boundary){
  # boundary operator for connection vectors: y-x
  # This changes the connection vector to respect periodic
  # boundary conditions on [0,1]x[0,1] so that the norm
  # of the connection vector is the shortest path in
  # periodic space
  boundary_operator_0_1 = function(x){
    ((x+.5) %% 1) - .5
  }
  dat = fuels$dat
  for(j in 1:fuels$reps){
    n = nrow(dat[[j]])
    XY = as.matrix(dat[[j]][,c('X','Y')])
    if(n>1){
      # D = as.matrix(dist(dat[[j]]))
      repulsion_mat = matrix(0,nrow=nrow(XY),ncol=2)
      for(i in 1:nrow(XY)){
        ref_point = as.numeric(XY[i,])
        
        # get connection vectors
        connect = sweep(XY,2,ref_point)
        # apply boundary operator
        if(boundary=='periodic'){
          connect = apply(connect,2,boundary_operator_0_1)
        }
        # get norm of connection vectors
        D = sqrt(rowSums(connect^2))
        # which points are B close to ref_point (excluding ref_point itself)
        close = which(D<B & D>0)
        # points = dat[[j]][close,c('X','Y')]
        
        if(length(close>0)){
          connect = connect[close,,drop=F]
          
          for(kk in 1:nrow(connect)){
            # d = ref_point - points[kk,]
            d = connect[kk,]
            if(sum(d!=0)){
              repulsion_mat[i,] = repulsion_mat[i,] + as.numeric(-d/(sum(d^2)))
            }
          }
        }
      }
      
      dat[[j]]$X = dat[[j]]$X + eps*repulsion_mat[,1]
      dat[[j]]$Y = dat[[j]]$Y + eps*repulsion_mat[,2]
      # bring out of bounds points within bounds
      limit = .01
      dat[[j]]$X = pmax(limit,dat[[j]]$X); dat[[j]]$X = pmin(dat[[j]]$X,1-limit)
      dat[[j]]$Y = pmax(limit,dat[[j]]$Y); dat[[j]]$Y = pmin(dat[[j]]$Y,1-limit)
      # native scale
      dat[[j]]$Xnat = dat[[j]]$X * fuels$dimX
      dat[[j]]$Ynat = dat[[j]]$Y * fuels$dimY
    }
  }
  fuels$dat = dat
  return(fuels)
}
# Generate 'reps' fuel maps with parameters theta over a domain of size [0,dimX]x[0,dimY]
gen_data = function(theta, dimX, dimY, heterogeneity.scale = 1, 
                    X.locs = NULL, X.vals = NULL, Beta = 0, 
                    reps = 1, GP.init.size=32, seed = NULL, 
                    I.transform='exp',
                    logis.scale=.217622, parallel=F, 
                    repulsion = F, repulsion.eps = 1e-3, repulsion.B = 1/2, repulsion.boundary = 'periodic')
{
  if(!is.null(seed))
    set.seed(seed)
  
  rho = theta[1]
  mu_r = theta[2]
  s2_r = theta[3]
  lambda = theta[4]
  dispersion = theta[5]
  if(length(theta)>5){
    repulsion.eps = theta[6]
    if(repulsion.eps>0){
      repulsion = T
    }
  }
  # if(length(theta)>5){
  #   mu_h = theta[6]
  #   s2_h = theta[7]
  # } else{
  #   mu_h = NULL
  # }
  mu_h = NULL
  
  # do everything on [0,1] and multiple to [dimX,dimY] later
  W = spatstat.geom::owin(c(0,1),c(0,1), mask=matrix(TRUE, GP.init.size,GP.init.size))
  W$dimX = dimX; W$dimY = dimY
  
  if(!is.null(X.vals)){
    # interpolate XB to cells defined by W
    XB = matrix(0,nrow=length(W$xcol),ncol=length(W$yrow))
    for(k in 1:length(Beta)){
      pred = expand.grid(x=W$xcol,y=W$yrow)
      XB = XB + Beta[k]*pracma::interp2(x=X.locs[[k]]$x,
                                        y=X.locs[[k]]$y,
                                        Z=X.vals[[k]],
                                        xp=pred$x,
                                        yp=pred$y)
    }
    XB = t(XB)
    XB = spatstat.geom::as.im(XB,W = W)
  } else if(Beta != 0){
    # No X's given, do constant mean at W coords
    XB = spatstat.geom::as.im(Beta,W = W)
  } else{
    XB = 0
  }
  
  suppressWarnings({
    # large lengthscales cause this:
    # Warning: _ out of _ terms (_%) in FFT calculation of matrix square root were negative, and were set to zero. Range: [-0.644, 708]
    gamma = spatstat.random::rGRFgauss(W = W, mu = XB, var = heterogeneity.scale, scale = rho,nsim = reps)
  })
  if(reps==1){gamma = list(gamma)}
  
  # sample number of points
  if(dispersion != 1){
    # mean parameterized CMP
    n_plus = mpcmp::rcomp(reps,mu = lambda*dimX*dimY, nu = dispersion)
  } else{
    # regular poisson
    n_plus = rpois(reps, lambda*dimX*dimY)
  }
  
  
  dat = vector(mode="list", length=reps)
  
  if(I.transform=='logistic'){
    gamma = spatstat.geom::solapply(gamma,plogis,location=0,scale=logis.scale)
  } else if(I.transform=='exp'){
    gamma = spatstat.geom::solapply(gamma,exp)
  } else if(I.transform=='relu'){
    gamma = spatstat.geom::solapply(gamma,relu)
  }

  #### sample pixels directly ####
  for(i in 1:reps){
    nn = n_plus[i]
    if(!is.finite(nn))
      stop(paste("Unable to generate Poisson process with a mean of",
                 nn, "points"))
    if(nn>0){
      dx = gamma[[i]]$xstep/2
      dy = gamma[[i]]$ystep/2
      df = as.data.frame(gamma[[i]])
      npix = nrow(df)
      lpix = df$value
      lpix[lpix<0] = 0
      ii = sample.int(npix, size=nn, replace=TRUE, prob=lpix)
      xx = df$x[ii] + runif(nn, -dx, dx)
      yy = df$y[ii] + runif(nn, -dy, dy)
      dat[[i]] = data.frame(X=xx,Y=yy,Xnat=xx*dimX,Ynat=yy*dimY)
      dat[[i]]$r = truncdist::rtrunc(nn,'norm',a=0,b=Inf,mu_r,sqrt(s2_r))
      if(!is.null(mu_h)){
        dat[[i]]$h = truncdist::rtrunc(nn,'norm',a=0,b=Inf,mu_h,s2_h)
      }
    } else{
      dat[[i]] = data.frame(X=numeric(),Y=numeric(),r=numeric())
      if(!is.null(mu_h)){
        dat[[i]]$h = numeric()
      }
    }
  }
  
  # return parameters as they may be needed to calculate prior
  fuels = list(dat=dat,W=W,dimX=dimX,dimY=dimY,reps=reps,
               theta=theta,
               heterogeneity.scale = heterogeneity.scale, 
               X.locs = X.locs, X.vals = X.vals, Beta = Beta, 
               GP.init.size=GP.init.size, I.transform=I.transform,seed = seed, logis.scale=logis.scale, parallel=parallel,
               gamma=gamma,nshrub=n_plus)
  class(fuels) = c('list','fuelsgen')
  
  if(repulsion){
    fuels = repulsion(fuels,repulsion.eps,repulsion.B,repulsion.boundary)
  }
  
  return(fuels)
}

# csv.filenames: vector of filenames for repeat observations
load_data = function(csv.filenames, dimX, dimY, heterogeneity.scale = 1, 
                    X.locs = NULL, X.vals = NULL, Beta = NULL, 
                    GP.init.size=32, logis.scale=.217622)
{
  reps = length(csv.filenames)
  dat = list(reps)
  for(i in 1:reps){
    # cols: X,Y,Xnat,Ynat,r,h
    dat[[i]] = data.frame(read.csv(csv.filenames[i]))
  }
  return(list(dat=dat,dimX=dimX,dimY=dimY,reps=reps,
              heterogeneity.scale = heterogeneity.scale, 
              X.locs = X.locs, X.vals = X.vals, Beta = Beta, 
              GP.init.size=GP.init.size, logis.scale=logis.scale))
}