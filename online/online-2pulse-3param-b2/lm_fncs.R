library(LaplacesDemon)
standardize_days = function(dys, shft, scl){
  (dys-shft)/scl
}


generate_one_lm = function(params, d2, days = 1:40,d1=15){

  d2.normalized = standardize_days(log(d2-d1), shft=2.5, scl=0.3)
  d22.normalized = standardize_days(log(d2-d1)^2, shft = 6, scl =2)
  d1_ind = (days>=d1)
  d2_ind = (days>=d2)
  p1.vec = generate_pd1_stacked(d2, max(days))[days]
  day.effects = invlogit(params[['mu']])-params[['beta.io']]*p1.vec
  d1.effect = params[['beta.rt1']]*d1_ind
  d2.coef = params[['beta.rt1']]*(1 + params[['beta.d']]*d2.normalized*params[['beta.io']] + params[['beta.d2']]*(d22.normalized)*params[['beta.io']])
  d2.effect = d2.coef*d2_ind
  #cat("\n",d2.coef, d2.normalized, d22.normalized)
  other.effects = d1.effect+d2.effect
  ltv = (day.effects)*days + other.effects
  ltv
}

generate_one_lm_inc = function(params, days = 1:20,d1=15){
  d1_ind = (days>=d1)
  p1.vec = generate_pd1_stacked(c(), max(days))[days]
  day.effects = params[['beta.io']]*p1.vec+ invlogit(params[['mu']])
  other.effects = params[['beta.rt1']]*d1_ind
  ltv = (day.effects)*(days) + other.effects
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

get_rteffects = function(one.pair){
  #d2 = min(which(one.pair[,'d']==2))
  days = one.pair[,'day']
  one.seq.fit = lm(one.pair[,'ltv']~days*one.pair[,'d']*one.pair[,'p']-1)
  one.seq.fit$coefficients
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
