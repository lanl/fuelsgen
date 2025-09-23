# This script contains user-facing functions should the fuels generation process be run
# from outside the shiny application.

#' @title Generate fuel layouts
#'
#' @description Generate fuel layouts on a rectangular domain using a handful of parameters to a generative model.
#' @param dimX (0,Inf): The size of the domain in X direction (meters)
#' @param dimY (0,Inf): The size of the domain in X direction (meters)
#' @param density (0,Inf): The mean number of fuels items
#' @param radius (0,Inf): The Average radius of a circular fuel element
#' @param sd_radius (0,Inf): The standard deviation of fuel element radii 
#' @param height (0,Inf): The Average height of a fuel element. Default is NULL indicating we should not sample heights.
#' @param sd_height (0,Inf): The standard deviation of fuel element heights. Default is 0 indicating no variance in heights.
#' @param heterogeneity (0,Inf): How heterogeneous the fuel maps will be. Larger values result in more heterogeneity.
#' @param heterogeneity.scale (0,Inf): How big of an effect should the heterogeneity have? Scaling factor.
#' @param reps (integer): The number of fuel layouts to generate
#' @param seed (integer): Optional random seed for reproducibility
#' @export
#'
gen_fuels = function(dimX, dimY,
                     density, dispersion = 1,
                     radius, sd_radius, 
                     height = NULL, sd_height = 0, 
                     heterogeneity, heterogeneity.scale = 1,
                     X.locs = NULL, X.vals = NULL, Beta = 0,
                     reps=1, GP.init.size = 32, seed=NULL, 
                     I.transform='exp',logis.scale=.217622,
                     repulsion=F, repulsion.eps = 1e-3, repulsion.B=1/2, repulsion.boundary='periodic',
                     parallel=F, verbose=T){
    if(verbose){
      cat('Generating',reps,'fuel maps.\n')
      if(is.null(seed))
        cat('For reproducible maps, set the random seed.\n')
      cat('Parameters:\n')
      cat('  density:       ',density,'\n')
      cat('  dispersion:    ',dispersion,'\n')
      cat('  fuel radius:   ',radius,'\n')
      cat('  fuel radius sd:',sd_radius,'\n')
      cat('  fuel height:   ',height,'\n')
      cat('  fuel height sd:',sd_height,'\n')
      cat('  heterogeneity: ',heterogeneity,'\n')
      cat('  heterogeneity scale: ',heterogeneity.scale,'\n')
      cat('  repulsion: ', repulsion,'\n')
      cat('  repulsion eps: ',repulsion.eps,'\n')
      cat('  repulsion B: ',repulsion.B,'\n')
    }
    
    theta = c(heterogeneity, radius, sd_radius^2, density/dimX/dimY, dispersion)
    if(repulsion){
      theta = c(theta,repulsion.eps)
    } else{
      theta = c(theta,0)
    }
    if(!is.null(height))
        theta = c(theta, height, sd_height)
    data = gen_data(theta, dimX, dimY, heterogeneity.scale, X.locs, X.vals, Beta, reps, GP.init.size, seed, I.transform, logis.scale, parallel, repulsion, repulsion.eps, repulsion.B, repulsion.boundary)
    return(data)
}

#' @title Plot fuel layouts
#'
#' @description Plot fuel layouts on a square tiled figure
#' @param data: fuels data object returned from gen_fuels
#' @export
#'
plot.fuelsgen = function(data,axis=F,circles=T,text.size=15,point_size=1,which=NULL,plot.dim=NULL,obs=NULL,fill='forestgreen'){
    if(data$reps>0){
      # make obs point patterns red
      cols = rep('black',data$reps)
      if(!is.null(obs)){
        cols[obs] = 'maroon'
      }
      if(!is.null(which)){
        data$dat = data$dat[which]
        data$reps = length(which)
      }
      if(is.null(plot.dim)){
        plot.dim = rep(min(10,ceiling(sqrt(data$reps))),2)
      }
      
      plot.list = vector(mode='list',length=data$reps)
      hmin = 0; hmax = 0

      for(i in 1:data$reps){
          if(!is.null(data$dat[[i]]$h)){
            tmp = max(data$dat[[i]]$h)
            if(tmp>hmax)
              hmax = tmp 
            tmp = min(data$dat[[i]]$h)
            if(tmp<hmin)
              hmin = tmp 
              plot.list[[i]] = ggplot2::ggplot() +
                  ggplot2::theme_bw() + 
                  ggplot2::coord_fixed(xlim = c(-.1, data$dimX+.1), ylim = c(-.1, data$dimY+.1))
              if(circles){
                plot.list[[i]] = plot.list[[i]] + ggforce::geom_circle(ggplot2::aes(x0 = Xnat, y0 = Ynat, r = r, fill = h, linewidth=I(.001)), data=data$dat[[i]])
              } else{
                plot.list[[i]] = plot.list[[i]] + ggplot2::geom_point(ggplot2::aes(x = Xnat, y = Ynat), size=point_size, data=data$dat[[i]], col=cols[i])
              }
              if(axis){
                plot.list[[i]] = plot.list[[i]] + ggplot2::theme(plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                               axis.text = ggplot2::element_text(size = text.size), 
                               axis.title = ggplot2::element_text(size = text.size),
                               axis.ticks = ggplot2::element_line(linewidth = 1.5))
              } else{
                plot.list[[i]] = plot.list[[i]] + ggplot2::labs(x='X (m)',y='Y (m)',fill='height (m)') + 
                                                  ggplot2::theme(plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                                                 axis.text = ggplot2::element_blank(), 
                                                                 axis.title = ggplot2::element_blank(),
                                                                 axis.ticks = ggplot2::element_blank())
              }
          } else{
              hmax = 0
              plot.list[[i]] = ggplot2::ggplot() +
                ggplot2::theme_bw() + 
                ggplot2::coord_fixed(xlim = c(-.1, data$dimX+.1), ylim = c(-.1, data$dimY+.1))
              if(circles){
                plot.list[[i]] = plot.list[[i]] + ggforce::geom_circle(ggplot2::aes(x0 = Xnat, y0 = Ynat, r = r,linewidth=I(.001)), data=data$dat[[i]],fill=fill)
              } else{
                plot.list[[i]] = plot.list[[i]] + ggplot2::geom_point(ggplot2::aes(x = Xnat, y = Ynat), size=point_size, data=data$dat[[i]], col=cols[i])
              }
              if(axis){
                plot.list[[i]] = plot.list[[i]] + ggplot2::theme(plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                                                 axis.text = ggplot2::element_text(size = text.size), 
                                                                 axis.title = ggplot2::element_text(size = text.size),
                                                                 axis.ticks = ggplot2::element_line(linewidth = 1.5))
              } else{
                plot.list[[i]] = plot.list[[i]] + ggplot2::labs(x='X (m)',y='Y (m)') + 
                  ggplot2::theme(plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                 axis.text = ggplot2::element_blank(), 
                                 axis.title = ggplot2::element_blank(),
                                 axis.ticks = ggplot2::element_blank())
              }
          }
      }
      
      if(!is.null(data$dat[[1]]$h)){
        return(patchwork::wrap_plots(plot.list,ncol=plot.dim[1],nrow=plot.dim[2], guides='collect') &
                 ggplot2::scale_fill_continuous(limits = c(hmin, hmax),breaks = seq(hmin,hmax,length.out=5),
                                                labels = function(x) sprintf("%.2f", x)))
      } else{
        return(patchwork::wrap_plots(plot.list,ncol=plot.dim[1],nrow=plot.dim[2], guides='collect'))
      }
  } else{
      # single map
      if(!is.null(data$dat$h)){
          myplot = ggplot2::ggplot() +
                   ggforce::geom_circle(ggplot2::aes(x0 = X, y0 = Y, r = r, fill = h), data=data$dat[[1]]) + 
                   ggplot2::labs(x='X (m)',y='Y (m)',fill='height (m)') + 
                   ggplot2::theme_bw() + 
                   ggplot2::theme(aspect.ratio = data$dimY/data$dimX) +
                   ggplot2::coord_fixed(xlim = c(-.1, data$dimX+.1), ylim = c(-.1, data$dimY+.1))
      } else{
          myplot = ggplot2::ggplot() + 
                   ggforce::geom_circle(ggplot2::aes(x0 = X, y0 = Y, r = r), fill=adjustcolor('grey',alpha.f=.5), data=data$dat[[1]]) + 
                   ggplot2::labs(x='X (m)',y='Y (m)') + 
                   ggplot2::theme_bw() + 
                   ggplot2::theme(aspect.ratio = data$dimY/data$dimX) +
                   ggplot2::coord_fixed(xlim = c(-.1, data$dimX+.1), ylim = c(-.1, data$dimY+.1))
      }
      # plot code for IF we move to ellipses instead of circles
      #data$dat$angle = runif(n=nrow(data$dat),min=0,max=360)
      #myplot = ggplot2::ggplot() + 
      #         ggforce::geom_ellipse(ggplot2::aes(x0 = X, y0 = Y, 
      #                                            a = .8*r, b=1.2*r, 
      #                                            angle=angle, fill=h), data=data$dat) + 
      #         ggplot2::labs(x='X (m)',y='Y (m)',fill='height (m)')
      return(myplot)
  }
}

#' @title Diagnostic plots for mcmc calibration
#'
#' @description Plot a fuelsgen_mcmc object
#' @param mcmc: fuelsgen_mcmc object from fuelsgen::mcmc_MH_adaptive()
#' @param burn: Number of initial samples to discard in plots
#' @param type: 'trace' or 'density'. 'trace' gives trace plots for the parameters, while 'density' plots the KDE and prior over the prior range
#' @export
#'
plot.fuelsgen_mcmc = function(mcmc,burn=0,type='trace',plot_prior=F,plot_llh=F,truth=NULL){
  par(mfrow=c(3,3)) # 6 parameters in the calibration
  names = mcmc$par_names
  names[names=='sigma']='radius variance'
  nsamp = nrow(mcmc$trace)
  lb = mcmc$prior$prior_params$lb
  ub = mcmc$prior$prior_params$ub
  # convert sigma to variance for prior scale
  # sig_id = which(names=='s2')
  # mcmc$trace[,sig_id] = 1/(mcmc$trace[,sig_id]^2)
  # mcmc$trace[,sig_id] = mcmc$trace[,sig_id]^2
  # ub[sig_id] = 1.1*max(mcmc$trace[,sig_id])
  # names[sig_id] = 'precision'
  if(type=='trace'){
    for(i in 1:length(names)){
      plot(mcmc$trace[(burn+1):nsamp,i],type='l',
           ylab=names[i])
      if(!is.null(truth)){abline(h=truth[i],lty=2,col='red')}
    }
    plot(mcmc$llh[(burn+1):nsamp],type='l',ylab='llh or D (ABC)')
  } else if(type=='density'){
    for(i in 1:length(names)){
      plot(density(mcmc$trace[(burn+1):nsamp,i],
                   from = lb[i],to = ub[i]),
           xlab=names[i],main='',ylab='density')
      if(!is.null(truth)){abline(v=truth[i],lty=2,col='red')}
      # abline(v=mcmc$prior$theta_est[i])
      if(plot_prior){
        switch(mcmc$prior$prior_params$dist[i],
               gamma = {
                 curve(dgamma(x,mcmc$prior$prior_params[[i]]$params[1],mcmc$prior$prior_params[[i]]$params[2]),add=T,lty=2)
               },
               uniform = {
                 curve(dunif(x,mcmc$prior$prior_params[[i]]$params[1],mcmc$prior$prior_params[[i]]$params[2]),add=T,lty=2)
               },
               truncnorm = {
                 curve(truncdist::dtrunc(x,'norm',mcmc$prior$prior_params[[i]]$params[3],mcmc$prior$prior_params[[i]]$params[4],
                                         mcmc$prior$prior_params[[i]]$params[1],mcmc$prior$prior_params[[i]]$params[2]),add=T,lty=2)
               },
               hcauchy = {
                 curve(LaplacesDemon::dhalfcauchy(x,mcmc$prior$prior_params[[i]]$params[1]),add=T,lty=2)
               })
      }
      # lines(priorfunc)
    }
    if(plot_llh)
      plot(density(mcmc$llh[(burn+1):nsamp]),xlab='llh or D (ABC)',ylab='density',main='')
  } else{
    stop('Currently supported types are trace and density')
  }
}

#' @title Plots for calibration priors
#'
#' @description Plot a fuelsgen_prior object
#' @param prior: fuelsgen_prior object from fuelsgen::get_prior_info()
#' @export
#'
plot.fuelsgen_prior = function(prior){
  par(mfrow=c(3,2)) # 5 parameters in the calibration
  lb = prior$lb
  ub = prior$ub
  names = prior$names
  for(i in 1:length(names)){
    switch(prior$dist[i],
           gamma = {
             curve(dgamma(x,prior$prior_params[[i]]$params[1],prior$prior_params[[i]]$params[2]),from=lb[i],to=ub[i],lty=2,ylab = 'density',xlab=names[i])
           },
           uniform = {
             curve(dunif(x,prior$prior_params[[i]]$params[1],prior$prior_params[[i]]$params[2]),from=lb[i],to=ub[i],lty=2,ylab = 'density',xlab=names[i])
           },
           truncnorm = {
             curve(truncdist::dtrunc(x,'norm',prior$prior_params[[i]]$params[3],prior$prior_params[[i]]$params[4],
                                     prior$prior_params[[i]]$params[1],prior$prior_params[[i]]$params[2]),
                   from=lb[i],to=ub[i],lty=2,ylab = 'density',xlab=names[i])
           },
           hcauchy = {
             curve(LaplacesDemon::dhalfcauchy(x,prior$prior_params[[i]]$params[1]),from=lb[i],to=ub[i],lty=2,ylab = 'density',xlab=ifelse(names[i]=='sigma','radius variance',names[i]))
           })
    abline(v=prior$theta_est[i],col='red',lty=2)
  }
}

#' @title Fix raster
#'
#' @description Change raster extent and coordinates to work with the generative model. Also scale raster values to lie in [-1,1]. Return a list containing the raster, coordinates, and values all in the correct format for the model. 
#' @param data: raster object
#' @export
#'
modify_raster = function(data,type='regular'){
  dimX = data@extent@xmax-data@extent@xmin
  dimY = data@extent@ymax-data@extent@ymin
  raster::extent(data) = c(0,dimX,0,dimY)
  resolution = raster::res(data)
  coords = raster::coordinates(data)
  locs = list()
  locs$x = unique(coords[,1])
  locs$y = rev(unique(coords[,2])) # y's must be in increasing order for pracma::interp2
  vals = raster::values(data,format='matrix')
  if(type=='regular'){
    vals[vals<0] = 0
    vals = (vals - min(vals)) / (max(vals) - min(vals))
    vals = 2*vals - 1 # scale to [-1,1] so that low canopy reduces probability of shrub
  }
  raster::values(data) = vals
  # flip y dim to match locs$y
  vals = vals[nrow(vals):1,]
  
  return(list(raster=data,dimX=dimX,dimY=dimY,locs=list(locs),vals=list(vals)))
}