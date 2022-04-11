setwd("~/RPAI/online/online-2pulse-working")
source("data_generation_fncs.R")
source("env_fncs.R")

set.seed(1998)
num_free_pulses=1
total_time = 40
wait_time=9
potential_actions = 1:11+wait_time

minibatch_size = 2
burn_in=500
inc.time = 15 - 2 + wait_time

all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
parameter_mat = make_parameter_mat(burn_in)
action_mat = matrix(sample(potential_actions, burn_in, replace=T),byrow=T, nrow=burn_in, ncol=num_free_pulses)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)

all_data= one_batch
all_actions = action_mat
all_parameters = parameter_mat
max_animals=700
animal.idx=nrow(all_data)
char.labs = get_labels(all.action.mat = all.action.mat)
q.fit=vector(mode="list", length(char.labs))
names(q.fit)=char.labs

while(animal.idx<max_animals){
  cat("animal:", animal.idx,"\n")
  parameter_mat = make_parameter_mat(minibatch_size)
  one_batch = generate_one_initial_batch(parameter_mat, inc.time, total_time)
  q.fit = update_model(q.fit,all_data, one_batch, action_mat)
  prediction.mat = get_predictions(q.fit, one_batch, all.action.mat)
  best_actions = names(q.fit)[apply(prediction.mat, c(1), which.max)]
  action_mat = as.matrix(as.numeric(best_actions), ncol=num_free_pulses)
  one_updated_batch = generate_one_batch(parameter_mat, action_mat, current.time = total_time)
  all_data = rbind(all_data, one_updated_batch)
  all_parameters = rbind(all_parameters, parameter_mat)
  all_actions = rbind(all_actions, action_mat)
  animal.idx = nrow(all_data)
}
