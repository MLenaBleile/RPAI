
generate_one_lm = function(params, d2, days = 1:40,d1=15){
  b2d2 = params[['rho']]+params[['bd2']]
  b2d1 = params[['rho']]+params[['bd1']]
  day.effects = (1 - d1*abs(params[['bd1']]) -(d1^2)*abs(params[['bd2']]) + d2*b2d1 +(d2^2)*b2d2)
  ltv = params[['b0']] + params[['bt']]*days*day.effects
  ltv
}

generate_one_lm_inc = function(params, days = 1:20,d1=15){
  day.effects = (1 - d1*abs(params[['bd1']]) -(d1^2)*abs(params[['bd2']]))
  ltv = params[['b0']] + params[['bt']]*days*day.effects
  ltv
}

get_loss =function(one.pair, params){
  d2 = min(which(one.pair[,'d']==2))
  days = one.pair[,'day']
  one.sim.sequence = generate_one_lm(params, d2=d2, days = days)
  one.sim.fit = lm(one.sim.sequence~days*one.pair[,'d']*one.pair[,'p'])
  one.seq.fit = lm(one.pair[,'ltv']~days*one.pair[,'d']*one.pair[,'p'])
  mean((one.seq.fit$coefficients-one.sim.fit$coefficients)^2)
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
