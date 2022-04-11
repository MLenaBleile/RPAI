setwd("~/RPAI/online/online-2pulse-3param-b1")
source("data_generation_fncs.R")
source("env_fncs.R")

set.seed(1998)
num_free_pulses=1
total_time = 40
wait_time=9
potential_actions = 1:11+wait_time
eps=1
eps_decay = .99
eps.vec=c(eps)
minibatch_size = 2
burn_in=50
inc.time = 15 - 2 + wait_time

all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
parameter_mat = make_parameter_mat(burn_in)
action_mat = matrix(sample(potential_actions, burn_in, replace=T),byrow=T, nrow=burn_in, ncol=num_free_pulses)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)

all_data= one_batch
all_actions = action_mat
all_parameters = parameter_mat
max_animals=1000
animal.idx=nrow(all_data)
char.labs = get_labels(all.action.mat = all.action.mat)
q.fit=vector(mode="list", length(char.labs))
names(q.fit)=char.labs


while(animal.idx<max_animals){
  cat("animal:", animal.idx,"\n")
  parameter_mat = make_parameter_mat(5)
  one.pair = generate_one(c(), parameter_mat[1,], current.time)
  q.fit = update_model(q.fit,all_data, one.pair, action_mat)
  prediction.mat = get_predictions(q.fit, one_batch, all.action.mat, add_variability = F)
  best_actions = names(q.fit)[apply(prediction.mat, c(1), which.min)]
  action_mat = as.matrix(as.numeric(best_actions), ncol=num_free_pulses)
  random.idx = which(runif(nrow(action_mat))>eps)
  action_mat[random.idx,] = sample(potential_actions, length(random.idx), replace=T)
  eps=eps*eps_decay
  eps.vec=c(eps.vec, eps)
  one_updated_batch = generate_one_batch(parameter_mat, action_mat, current.time = total_time)
  all_data = rbind(all_data, one_updated_batch)
  all_parameters = rbind(all_parameters, parameter_mat)
  all_actions = rbind(all_actions, action_mat)
  animal.idx = nrow(all_data)
}
