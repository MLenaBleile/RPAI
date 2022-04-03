setwd("~/Documents/Dissertation/RPAI/online/online-3pulse-3param-v7")
source("data_generation_fncs.R")
source("env_fncs.R")

set.seed(1998)
max_epochs=5000
test_num=500
parameter_mat = make_parameter_mat(max_epochs+test_num)

total_time=60
num_pulses=3
num_free_pulses=num_pulses-1
state_size = total_time+11+num_pulses
#state_size = total_time+num_pulses
state_size=9+num_pulses

waitime_vec=rep(8, num_free_pulses)
state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
potential_actions = 1:12
action_size=length(potential_actions)

q.fit=readRDS("model.rds")

test_indices=1:500
num_mice = length(test_indices)
optimal_actions=rep(NA, num_mice)
data_mat = matrix(NA, ncol=20, nrow=num_mice)
true.mins = rep(NA, num_mice)
agent.outcomes=rep(NA, num_mice)
day15.outcomes=rep(NA, num_mice)
selected_actions=c()

agent.action.mat = matrix(NA, nrow=num_mice, ncol=num_free_pulses)
reference_days=c(8)

agent.outcome.vec=rep(NA, 500)
refdays = c(10,14,20)
ref.outcome.mat = matrix(NA, nrow=num_mice, ncol=length(refdays)+1)
references=c(paste("day", refdays), "random")
colnames(ref.outcome.mat) = c(paste("day", refdays), "random")

for(mouse in 1:500){
  parameter_vec=parameter_mat[mouse,]
  action.vec=c()
  one.sequence = generate_one(action.vec, parameter_vec, maxtime=20)
  while(length(action.vec)<num_free_pulses){
    isdone=(length(action.vec)==(num_free_pulses-1))
    one.state=sequence_to_state(one.sequence = as.numeric(one.sequence), action.vec = action.vec, done=isdone, num_free_pulses = 2)
    one.action = get_action(q.fit,one.state, potential_actions = 1:10)+waitime_vec[length(action.vec)+1]
    action.vec=c(action.vec, one.action)
    seqlength = (20+sum(action.vec))*(1-isdone)+60*isdone
    one.sequence = generate_one(cumsum(action.vec)+15, parameter_vec, maxtime = seqlength)
  }
  for(ref.idx in 1:length(refdays)){
    ref = refdays[ref.idx]
    one.reference = generate_one(cumsum(c(ref,ref))+15, parameter_vec, maxtime = seqlength)
    ref.outcome.mat[mouse,ref.idx]=log(one.reference[seqlength])
  }
  one.random = sample(1:10+waitime_vec[length(action.vec)],2, replace=T)
  ref.outcome.mat[mouse,length(refdays)+1]= log(generate_one(cumsum(one.random)+15, parameter_vec, maxtime = seqlength))[seqlength]
  agent.action.mat[mouse,]=action.vec
  agent.outcome.vec[mouse] = log(one.sequence[seqlength])
  
}
effect.sizes = rep(NA, ncol(ref.outcome.mat))
names(effect.sizes) = references

adj.val=2
plot(density(-agent.outcome.vec+ref.outcome.mat[,1], adjust=adj.val), col="red", xlab="final ltv reference-final ltv agent",main="3 pulse performance", xlim=c(-.7,.7), ylim=c(0,15))
colors= c("red","blue", "orange", "purple")
for(ref.idx in 1:ncol(ref.outcome.mat)){
  ref = refdays[ref.idx]
  lines(density(ref.outcome.mat[,ref.idx] - agent.outcome.vec, adjust=adj.val), col=colors[ref.idx])
  print(t.test(ref.outcome.mat[,ref.idx] - agent.outcome.vec))
  delta = ref.outcome.mat[,ref.idx] - agent.outcome.vec
  delta2 = -ref.outcome.mat[,ref.idx] + ref.outcome.mat[,'random']
  effect.sizes[[references[ref.idx]]] = mean(delta)/sd(delta)
  print(summary(delta))
  print(mean(delta[delta<=0]))
  print(mean(delta[delta>0]))
  cat("proportion of times agent was better than", colnames(ref.outcome.mat)[ref.idx],"day spacing:" )
  print(1-length(delta[delta<=0])/100)
  cat("proportion of times ", colnames(ref.outcome.mat)[ref.idx],"was better than random" )
  print(1-length(delta[delta2<=0])/100)
}
legend("topright", legend=paste(c(paste("day", refdays),"random"), round(effect.sizes, digits=3), sep=": d="), pch=16, col=colors)
abline(v=0, lty=2)
mean(-agent.outcome.vec+ref.outcome.mat[,1])
