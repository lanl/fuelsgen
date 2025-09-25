library(fuelsgen)
lambda.true = 100
rho.true = .1
mu.true = 1
sigma.true = .1
nu.true = 1
fuel = gen_fuels(dimX = 30, dimY = 30,
                 density = lambda.true,        #expected shrubs
                 dispersion = nu.true,
                 heterogeneity = rho.true,     # level of heterogeneity in shrub placement
                 heterogeneity.scale = 1,      # scale of the mean-zero GP realization
                 radius = mu.true, sd_radius = sigma.true,  # normal distribution parameters for shrub radius
                 height = NULL, sd_height = 0, # normal distribution parameters for shrub height
                 reps=10,                      # number of random maps to generate
                 GP.init.size=32,             # How densely to sample GP (n^3 scaling, <=1000 is pretty fast)
                 seed=10,                      # random seed for reproducibility
                 parallel=F)                   # Parallel option will be faster for expensive generation only (parallel overhead), rng seed doesn't work with parallel
plot(fuel)
info = get_mets_info(En = T, Vn = T)
y_obs = get_mets(fuel,info)
prior = get_prior_info(fuel,y_obs,est_cov_obs = T,est_cov_samples = 25, est_cov_reps = 25)
MCMC = mcmc_MH_adaptive(y_obs,fuel,prior,
                        n.samples=1000,n.burn=0,gen_reps=24,
                        adapt.par = c(100,20,.5,.75),
                        prop.sigma = diag(.2^2*prior$theta_est), 
                        ABC = T, ABC_eps = 12, ABC_wt = c(1,.5))
plot(MCMC)
scatter_pairs(list(MCMC$trace[501:1000,]),prior = prior, truth = prior$theta_est)
