setwd("C:/Users/s198663/Documents/RPAI/online/online-2pulse-working")
source("data_generation_fncs.R")
source("env_fncs_2pulse.R")

set.seed(1998)
max_animals=1000
test_num=100

maxtime = 40
num_free_pulses=1
eps=1
eps_decay=.99
eps.vec=c()

wait_time=7
state_mat = matrix(NA, nrow = 0, ncol=4)
colnames(state_mat) = c("ltv", "day", 'd','animal')
potential_actions = 1:13 + wait_time


minibatch_size=100



rewards = c()
mouse=0
epoch=1

potential.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
colnames(potential.action.mat)="action"

all_data_mat = matrix(nrow=0, ncol=4)
colnames(all_data_mat) = c("ltv", "day","animal", "action")

minibatch_parameters = make_parameter_mat(minibatch_size)
parameter_mat = minibatch_parameters

parameter_vec= minibatch_parameters[1,]
one.sequence=generate_one(c(), parameter_vec = parameter_vec, maxtime = 20)$ltv

data_mat = generate_one_batch(minibatch_size = minibatch_size, minibatch_parameters = minibatch_parameters, wait_time=wait_time)

random_action_idx = sample(1:nrow(potential.action.mat), minibatch_size, replace=T)
action_mat = matrix(potential.action.mat[random_action_idx,], ncol=num_free_pulses)
data_mat = update_one_batch(data_mat,minibatch_parameters, action_mat, maxtime=40)
all_data_mat = data_mat
animal.idx=max(all_data_mat[,'animal'])

while(animal.idx < max_animals){
  cat("animal idx: ", animal.idx,"\n")
  minibatch_parameters = make_parameter_mat(minibatch_size)
  parameter_mat = cbind(parameter_mat, minibatch_parameters)
  
  data_mat = generate_one_batch(minibatch_size = minibatch_size, minibatch_parameters = minibatch_parameters, wait_time=wait_time, animal.idx = animal.idx)
  data_for_model = prep_data_for_model(all_data_mat, data_mat, potential_actions = potential_actions)
  q.fit = update_model(data_for_model)
  eps =eps*eps_decay
  eps.vec = c(eps.vec, eps)
  action_mat = matrix(NA, nrow = minibatch_size, ncol=num_free_pulses)
  for(one.animal.idx in 1:minibatch_size){
    one.animal = one.animal.idx + animal.idx - minibatch_size
    if(rnorm(1)<eps & one.animal>minibatch_size){
    action_mat[one.animal.idx,] = get_action(q.fit=q.fit, one.animal = one.animal, maxtime=40, potential.action.mat = potential.action.mat, nsamples=1)
    }else{
      action_mat[one.animal.idx,] = sample(potential_actions,1)
    }
  }
  new_data_mat = update_one_batch(data_mat,minibatch_parameters = minibatch_parameters, action_mat = action_mat, maxtime=40)
  all_data_mat = rbind(all_data_mat, new_data_mat)
  animal.idx = max(all_data_mat[,'animal'])
}

