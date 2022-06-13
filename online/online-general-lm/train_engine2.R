setwd("~/RPAI/online/online-general-lm")
source("data_generation_fncs.R")
source("env_fncs.R")
source("lm_fncs.R")
library(glmnet)
library(lme4)

set.seed(1998)
num_free_pulses=1
total_time = 40 + 20*(num_free_pulses-1)
gen.mod = "sumexp"
wait_time=7

potential_actions = 1:5+wait_time
train_num=50/num_free_pulses
adapt_num=150/num_free_pulses
test_num = 50
train_idx=1:(train_num*num_free_pulses)
adapt_idx=1:(adapt_num*num_free_pulses) + max(train_idx)

test_idx = tail(c(adapt_idx),test_num*num_free_pulses)
total_mice=train_num+adapt_num
num_pc = 3
inc.time=wait_time+15
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
optim.acts=numeric(total_mice)

##make initial data
parameter_mat = make_parameter_mat(total_mice, gen.mod=gen.mod)
action_mat = all.action.mat[sample(1:nrow(all.action.mat), total_mice, replace=T),]
if(num_free_pulses==1){
  action_mat = matrix(action_mat, ncol=1)
}
action_mat[1:nrow(all.action.mat),] = as.matrix(all.action.mat)
one.batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time, gen.mod=gen.mod)
one.batch.molten = reshape2::melt(one.batch[,,c("ltv","d")], varnames=c("animal", "day","feature"))
one.batch.clean = dcast(one.batch.molten, animal+day~feature)
one.batch.clean$sqrtday = sqrt(one.batch.clean$day)
fit = lmer(ltv~ sqrtday + d + sqrtday:d|animal, data=one.batch.clean[1:train_num,])


for(one.mouse in 1:train_num){
  op.at=get.optim.plan(parameter_mat[one.mouse,],total_time, as.matrix(all.action.mat), gen.mod=gen.mod)
  optim.acts[one.mouse] = op.at
}
library(reshape2)

for(one.mouse in adapt_idx){
  parameter_vec = parameter_mat[one.mouse,]
  new.data = melt(generate_one_sumexp(c(), parameter_vec, inc.time)[,,c("ltv","d")], varnames=c("day","feature"))
  new.data$animal = one.mouse
  new.data$sqrtday = sqrt(new.data$day)
  predictions = numeric(length(potential_actions))
  
}
