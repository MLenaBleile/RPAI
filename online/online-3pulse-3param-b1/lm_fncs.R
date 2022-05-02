library(LaplacesDemon)
standardize_days = function(dys, shft, scl){
  (dys-shft)/scl
}



get_loss =function(one.pair, params){
  d2 = min(which(one.pair[,'d']==2))
  days = one.pair[,'day']
  one.sim.sequence = generate_one_lm(params, d2=d2, days = days)
  one.sim.fit = lm(one.sim.sequence~days*one.pair[,'d']*one.pair[,'p'])
  one.seq.fit = lm(one.pair[,'ltv']~days*one.pair[,'d']*one.pair[,'p'])
  mean((one.seq.fit$coefficients-one.sim.fit$coefficients)^2)
}

get_rteffects = function(one.pair){
  #d2 = min(which(one.pair[,'d']==2))
  days = one.pair[,'day']
  dose_norm = logit((one.pair[,'d']+.5)/4)
  one.seq.fit = lm(one.pair[,'ltv']~days*dose_norm*one.pair[,'p']-1)
  ss=summary(one.seq.fit)
  day.rescale = invlogit(max(one.pair[,'day'])/60)

  c(one.seq.fit$coefficients, summary(one.seq.fit)$r.squared)
}

get_loss_inc = function(params, one.pair){
  days = one.pair[,'day']
  one.sim.sequence = generate_one_lm_inc(params,days = days)
  one.sim.fit = lm(one.sim.sequence~days*one.pair[,'d'] + one.pair[,'d']:one.pair[,'p'])
  one.seq.fit = lm(one.pair[,'ltv']~days*one.pair[,'d'] + one.pair[,'d']:one.pair[,'p'])
  mean((one.seq.fit$coefficients-one.sim.fit$coefficients)^2)
}

obj.fnc=function(params, all_data){
  loss.vec=rep(NA, dim(all_data)[1])
  for(midx in 1:length(loss.vec)){
    one.pair=all_data[midx,,]
    loss.vec[midx] = get_loss(one.pair, params)
  }
  log(mean(loss.vec))
}
