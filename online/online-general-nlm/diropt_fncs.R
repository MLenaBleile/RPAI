##function get_fitted_params
##in: one.pair
##out: best fitting params
initialize_params = function(param_names){
  num_params=length(param_names)
  random.init = runif(1:num_params, para_set[param_names,'lower'], para_set[param_names,'upper'])
  names(random.init) = param_names
  random.init
}

r.to.unit = Vectorize(function(x, s=1){s/(1+exp(x))}, vectorize.args=c("x","s"))
unit.to.r = Vectorize(function(y,s=1){log(s/y - 1)}, vectorize.args = c("y","s"))

get_loss = function(params.transformed,one_batch, num_free_pulses, total_time){
  num_in_batch=dim(one_batch)[1]
  loss.vec = numeric(dim(one_batch)[1])
  para_set = make_default_para_set()
  params=r.to.unit(params.transformed, s=para_set[names(params.transformed),'upper'])
  for(one in 1:num_in_batch){
    one.pair = one_batch[one,,]
    
    observed.effects= get_rteffects(one.pair, num_free_pulses = num_free_pulses,total_time=total_time)
    rt.gy=one.pair[,'d']
    rt.idx = which(rt.gy[-1] != rt.gy[-length(rt.gy)])
    radiation_days = c(one.pair[rt.idx,'day']+1)
    simulated.pair = generate_one(radiation_days=radiation_days[-1], parameter_vec=params, maxtime=total_time)
    simulated.effects = get_rteffects(simulated.pair[1,1:max(one.pair[,'day']),], num_free_pulses = num_free_pulses, total_time=total_time)
    msq.e=mean((simulated.pair[1,1:max(one.pair[,'day']),'ltv']-one.pair[,'ltv'])^2)
    loss.vec[one] = sqrt(mean((simulated.effects-observed.effects)^2)) + sqrt(msq.e)
  }
  mean(loss.vec)
}

get_optim_params = function(one_batch, param_names, num_free_pulses, total_time, init.par=NA){
para_set=make_default_para_set()
if(is.na(init.par[1])){
init.par=initialize_params(param_names)}else{}
random.init=unit.to.r(init.par,s=para_set[param_names,'upper'])
optim.obj = optim(par = random.init,fn=get_loss,method="Nelder-Mead",one_batch=one_batch, num_free_pulses=num_free_pulses, total_time=total_time)
r.to.unit(optim.obj$par, s=para_set[param_names,'upper'])
}
