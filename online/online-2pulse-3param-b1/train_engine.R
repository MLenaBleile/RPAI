setwd("~/RPAI/online/online-2pulse-3param-b1")
source("data_generation_fncs.R")
source("env_fncs.R")
source("lm_fncs.R")

set.seed(1998)
num_free_pulses=1
total_time = 40
wait_time=9
potential_actions = 1:5+wait_time
eps=1
eps_decay = .995
eps.vec=c(eps)
minibatch_size = 2
burn_in=10
inc.time = 15 - 2 + wait_time

all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
parameter_mat = make_parameter_mat(burn_in)
action_mat = matrix(sample(potential_actions, burn_in, replace=T),byrow=T, nrow=burn_in, ncol=num_free_pulses)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)

all_data= one_batch
all_actions = action_mat
all_parameters = parameter_mat
max_animals=60
max_updates = 40
animal.idx=nrow(all_data)
char.labs = get_labels(all.action.mat = all.action.mat)
params=rep(0.5,5)
names(params) = c("b0","bt","bd1","bd2","rho")
dimnames(all_data) = list(1:animal.idx, 1:total_time, c("ltv","d","p","day"))


while(animal.idx<max_animals){
  cat("animal:", animal.idx,"\n")
  parameter_mat = make_parameter_mat(5)
  one.pair = generate_one(c(), parameter_mat[1,], inc.time)
  
  num.extra = -dim(one.pair)[2] + dim(all_data)[2]
  one.pair.extended = abind::abind(one.pair, array(NA, dim=c(1, num.extra,4)), along=2)
  dimnames(one.pair.extended) = list(animal.idx+1, 1:total_time,c("ltv","d","p","day"))
  if(animal.idx<max_updates){
  params = update_model(all_data)}else{}
  one.pair=one.pair[1,,]
  prediction.vec = get_predictions(params,one.pair, all.action.mat)
  best_action = all.action.mat[which.min(prediction.vec),]
  action_mat = as.matrix(as.numeric(best_action), ncol=num_free_pulses)
  #random.idx = which(runif(nrow(action_mat))<eps)
  parameter_mat=t(as.matrix(parameter_mat[1,], nrow=1, ncol=8))
  one_updated_batch = generate_one_batch(parameter_mat, action_mat, current.time = total_time)
  all_data = abind::abind(all_data, one_updated_batch,along=1)
  
  all_parameters = rbind(all_parameters, parameter_mat)
  all_actions = rbind(all_actions, action_mat)
  animal.idx = nrow(all_data)
  dimnames(all_data) = list(1:animal.idx, 1:total_time, c("ltv","d","p","day"))
}


