setwd("~/RPAI/online/online-2pulse-nn-recursive")
source("datagen_fncs.R")
source("env_fncs.R")
para_set=make_default_para_set()

set.seed(1998)
calibration_num=500
adapt_num = 500
total_num=calibration_num+adapt_num
gen.mod="recursive"

maxtime = 40
num_free_pulses=1
wait_time=7
potential_actions = 1:7 + wait_time
##weights for the raw ssq and effect ssq in loss
loss.weights=c(0,1,0)
names(loss.weights) = c("effects",'raw','liklihood')
overall.loss.weights = c(0,1,0)
names(overall.loss.weights) = c("effects",'raw','liklihood')
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
##this is the date of the first IO treatment (two days before the first RT), plus the waiting time.
inc_days = 15-2+wait_time
true_parameters=make_parameter_mat(total_num, gen.mod)

initial.sequences = matrix(nrow=total_num, ncol=inc_days)
all.pairs = array(dim=c(total_num,maxtime,4))
dimnames(all.pairs) = list(NULL, NULL, c("ltv","d","p",'day'))
names(dimnames(all.pairs)) = c('animal','time', 'feature')

#initialize a bunch of objects for storing the outcomes
agent.outcomes=rep(NA, total_num)
reference_days=c(1,potential_actions,20)
reference.outcomes=matrix(NA, nrow=total_num, ncol=length(reference_days)+1)
colnames(reference.outcomes) = c(paste("day",reference_days),'random')
selected_actions=matrix(NA, nrow=total_num, ncol=num_free_pulses)
optimal_actions=numeric(total_num)

mouse=1
for(mouse in mouse:total_num){
  cat("optimizing mouse ",mouse,"\n")
  true.param.vec= true_parameters[mouse,]
  one.reference.pairs = generate_one_counterfactualset(parameter_vec=true.param.vec,total_days = maxtime,potential_actions = 1:20, gen.mod=gen.mod)
  one.optima= which.min(one.reference.pairs[potential_actions,maxtime,'ltv'])
  optimal_actions[mouse] = one.optima+wait_time
  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
    #day15.outcomes[mouse] = one.reference[reference_day-6,total_days]
    reference.outcomes[mouse, refday.idx]=one.reference.pairs[reference_day,maxtime,'ltv']
  }
  
  random.action = sample(potential_actions, size=num_free_pulses, replace=T)
  reference.outcomes[mouse,'random'] = one.reference.pairs[random.action,maxtime,'ltv']
  inc.pair = generate_one(radiation_days = c(), parameter_vec = true.param.vec, maxtime=inc_days, gen.mod=gen.mod)
  initial.sequences[mouse,] = inc.pair[1,,'ltv']
  
  if(mouse<calibration_num){
    one.action=random.action
  }else{
    q.nn = update_q(q.nn,all.pairs, maxtime=maxtime, inc_days=inc_days)
    one.action = get.optim.plan(inc.pair=inc.pair,q.nn = q.nn, all.action.mat = all.action.mat)
    selected_actions[mouse,] = one.action
  }
  selected_actions[mouse,] = one.action
  one.agent.treated.sequence = generate_one(cumsum(one.action)+15, parameter_vec=true.param.vec, maxtime=maxtime, gen.mod=gen.mod)
  all.pairs[mouse,,] = one.agent.treated.sequence[1,,]
  agent.outcomes[mouse] = one.agent.treated.sequence[1,maxtime,'ltv']
  
}
