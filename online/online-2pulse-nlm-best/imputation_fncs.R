impute = function( one.inc.sequence,one.action.vec,samples, imptimes, impnum=1){
  known.idx = 1:(length(one.inc.sequence)+length(one.action.vec))
    #get a dataframe with just the least squares means
    midx=which(substr(colnames(samples), 1,3)=="lsm")
    meanss = samples[,midx]
    if(length(imptimes>0)){
      #create diagonal and off-diagonal matrix elements of the partitioned Sigma matrix (corresponding to known vs unknown obs)
      sig.offdiag = matrix(data=rep(samples$covp1[impnum],length(imptimes)*(ncol(meanss)-length(imptimes))), ncol=length(imptimes))
      p1 = ncol(meanss) - length(imptimes)
      p2=length(imptimes)
      sig.11 = matrix(data=rep(samples$covp1[impnum], p1*p1), ncol=p1)
      diag(sig.11) = rep(samples$covp2[impnum]+samples$covp1[impnum], p1)
      sig.22 = matrix(data=rep(samples$covp1[impnum], p2*p2), ncol=p2)
      diag(sig.22) = rep(samples$covp2[impnum]+samples$covp1[impnum], p2)
      #get the desired parameters for the distribtution
      muvec = meanss[impnum,imptimes] +t(sig.offdiag)%*%solve(sig.11)%*%(as.matrix(as.numeric(meanss[impnum, -imptimes])-as.numeric(one.inc.sequence)))
      
      sigmat = sig.22 - t(sig.offdiag)%*%solve(sig.11)%*%sig.offdiag
      #generate obs from this distribution and impute
      data[aidx,(imptimes+1)] = MASS::mvrnorm(n=1, mu=as.numeric(muvec), Sigma=sigmat)
      
    }
  
  data
}


###Function to generate samples using a mixed model
get_samples = function(one.group,days,prior_obj,idcols=1, nsamp=5, chains=4, cores=4,iter=1000,warmup=500){
  #names(one.group)[animalcol]="Animal"
  data.long = reshape(one.group,
                      direction="long",
                      varying = list(names(one.group)[-c(idcols)]),
                      v.names="LTV",
                      idvar= list(names(one.group)[idcols]),
                      times = days)
  tgmix = data.long
  #imp = mice::mice(tgmix)
  ##fit using time as a categorical like before
  tgmix_fit = brms::brm(data = tgmix,
                        family = gaussian,
                        prior = prior_obj,
                        formula = LTV ~ as.factor(days) + (1|animal),
                        iter = iter, warmup = warmup, chains = chains, cores = cores,
                        control = list(adapt_delta = .975, max_treedepth = 20),
                        seed = 190831)
  samples_init = brms::posterior_samples(tgmix_fit, subset = sample(warmup:iter, nsamp))
  get.mus1 = function(colnum){samples_init$b_Intercept+ samples_init[,colnum]}
  post.mus = sapply(2:length(days), get.mus1)
  post.mus = cbind(as.matrix(samples_init[,1]), post.mus)
  colnames(post.mus) = paste("lsm",1:length(days),sep="")
  
  covpTime = (samples_init$sigma)^2
  
  
  covpAnimal = (samples_init$sd_Animal__Intercept)^2
  
  samples = cbind(as.matrix(covpTime),as.matrix(covpAnimal),post.mus)
  colnames(samples)[1:2] = c("covp2","covp1")
  
  samples = as.data.frame(samples)
  mono.inc.idx= apply(samples[,-c(1,2)], 1, function(x)all(x==cummax(x)))
  samples = samples[all()]
  samples$sample = 1:nsamp
  samples
}