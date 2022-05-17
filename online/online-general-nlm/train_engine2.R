setwd("~/RPAI/online/online-general-nlm")
source("data_generation_fncs.R")
source("env_fncs.R")
source("lm_fncs.R")
source("diropt_fncs.R")
para_set=make_default_para_set()

set.seed(1998)
num_free_pulses=1
total_time = 40 + 20*(num_free_pulses-1)
wait_time=9
potential_actions = 1:6+wait_time
train_num=30/num_free_pulses
adapt_num=120/num_free_pulses
test_num = 70
train_idx=1:(train_num*num_free_pulses)
adapt_idx=1:(adapt_num*num_free_pulses) + max(train_idx)
nonadapt_idx=1:(test_num*num_free_pulses) + max(adapt_idx)
test_idx = tail(c(adapt_idx),test_num*num_free_pulses)
total_mice=train_num+adapt_num
inc.time=wait_time+15-2
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)

##make initial data
parameter_mat = make_parameter_mat(total_mice)
action_mat = all.action.mat[sample(1:nrow(all.action.mat), total_mice, replace=T),]
if(num_free_pulses==1){
  action_mat = matrix(action_mat, ncol=1)
}
action_mat[1:nrow(all.action.mat),] = as.matrix(all.action.mat)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)

one.pair= one_batch[1,1:total_time,]
test_effects = get_rteffects(one.pair, num_free_pulses, total_time)
rteffects = array(0, dim=c(total_mice, length(test_effects), num_free_pulses))
dimnames(rteffects) = list(1:total_mice,c(names(test_effects)), as.character(1:num_free_pulses+1))
names(dimnames(rteffects)) = c("animal", "effect", "pulse")

##calibrate on initial data and record true optimal actions
optim.acts= matrix(nrow=total_mice, ncol=num_free_pulses)

dn.mat = matrix(NA, nrow=total_mice, ncol=num_free_pulses)
param_names = c('mu','rho','lambda')
estimated.parameter.mat = matrix(nrow=adapt_num, ncol =length(param_names))
colnames(estimated.parameter.mat) = param_names
for(one.mouse in 1:train_num){
  op.at=get.optim.plan(parameter_mat[one.mouse,],total_time, as.matrix(all.action.mat))
  optim.acts[one.mouse,] = op.at
  
  one.pair=one_batch[one.mouse,,]
  
  one.effects = get_rteffects(one_batch[one.mouse,1:inc.time,],num_free_pulses = num_free_pulses, total_time=inc.time)
  rteffects[one.mouse,,1] = one.effects
  #one.paraset = get_optim_params(one.pair,param_names = param_names, num_free_pulses=num_free_pulses, total_time=total_time)
  #estimated.parameter.mat[one.mouse,param_names] = one.paraset[param_names]
}

est.par.means=colMeans(estimated.parameter.mat, na.rm=T)

predicted.optims=matrix(nrow=adapt_num, ncol=num_free_pulses)
overall.means = get_optim_params(one_batch = one_batch[1:one.mouse,,],param_names = param_names, num_free_pulses=num_free_pulses, total_time=total_time)
###do adaptive part
fixed.params=c()
weight.factors=matrix(NA, nrow=adapt_num, ncol=ncol(estimated.parameter.mat))
for(one.mouse in train_num+1:adapt_num){
  cat("\n mouse:", one.mouse,"\n")
  op.at=get.optim.plan(parameter_mat[one.mouse,],total_time, as.matrix(all.action.mat))
  optim.acts[one.mouse,] = op.at
  one.pair=generate_one(c(), parameter_mat[one.mouse,], maxtime=inc.time)
  one.paraset = estimated.parameter.mat[one.mouse-train_num,]
    #get_optim_params(one.pair,param_names = param_names, num_free_pulses=num_free_pulses, total_time=total_time)
  one.effects = get_rteffects(one_batch[one.mouse,1:inc.time,],num_free_pulses = num_free_pulses, total_time=inc.time)
  rteffects[one.mouse,,1] = one.effects
  
  weight.factor.raw = sum(abs(one.pair[1,1:inc.time,'ltv']-colMeans(one_batch[1:one.mouse,1:inc.time,'ltv'])))
  weight.factor=min(sqrt(weight.factor.raw)*2,1)
  weight.factors[one.mouse-train_num,]=weight.factor
  #weight.factor=.9
  
  #if((one.mouse-1)%%10==0){
  #overall.means = get_optim_params(one_batch = one_batch[1:one.mouse,,],param_names = param_names, init.par=overall.means,num_free_pulses=num_free_pulses, total_time=total_time)}else{}
  predictive.params = one.paraset*(weight.factor)+overall.means*(1-weight.factor)
  
  if(max(weight.factor)>1){
    cat("\nalert: weight factor=", weight.factor,"\n")
    beepr::beep()
  }
  predictive.params[fixed.params] = overall.means[fixed.params]
  predicted.optim = get.optim.plan(predictive.params, total_time, as.matrix(all.action.mat))
  predicted.optims[one.mouse-train_num,]=predicted.optim
  one.pair = generate_one(c(cumsum(predicted.optim)+15), parameter_mat[one.mouse,], maxtime=total_time)
  one_batch[one.mouse,,] = one.pair[1,,]
  #one.paraset = get_optim_params(one.pair,param_names = param_names, num_free_pulses=num_free_pulses, total_time=total_time)
  estimated.parameter.mat[one.mouse-train_num,param_names] = one.paraset[param_names]
  
}
beepr::beep()



source("test_engine.R")
