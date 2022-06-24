library(LaplacesDemon)
standardize_days = function(dys, shft, scl){
  (dys-shft)/scl
}



get_loss =function(one.pair, params){
  d2 = min(which(one.pair[,'d']==2))
  days = one.pair[,'day']
  one.sim.sequence = generate_one_lm(params, d2=d2, days = days)
  one.sim.fit = lm(one.sim.sequence~days*one.pair[,'d'])
  one.seq.fit = lm(one.pair[,'ltv']~days*one.pair[,'d'])
  mean((one.seq.fit$coefficients-one.sim.fit$coefficients)^2)
}

get_rteffects = function(one.pair, num_free_pulses, total_time){
  #d2 = min(which(one.pair[,'d']==2))
  print(one.pair[,'ltv'])
  days = sqrt(one.pair[,'day'])
  dose_norm = logit((one.pair[,'d']+.5)/(num_free_pulses+2))
  one.seq.fit = lm(one.pair[,'ltv']~days*dose_norm)
  ss=summary(one.seq.fit)
  #day.rescale = invlogit(max(one.pair[,'day'])/total_time)

  out.vec=c(one.seq.fit$coefficients, summary(one.seq.fit)$r.squared)
  names(out.vec)[length(out.vec)] = 'rsq'
  out.vec
}

get_loss_inc = function(params, one.pair){
  days = one.pair[,'day']
  one.sim.sequence = generate_one_lm_inc(params,days = days)
  one.sim.fit = lm(one.sim.sequence~days*one.pair[,'d'] )
  one.seq.fit = lm(one.pair[,'ltv']~days*one.pair[,'d'] )
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
