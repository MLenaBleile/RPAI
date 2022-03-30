setwd("~/Documents/Dissertation/RPAI/online/online-4pulse-2param")
source("data_generation_fncs.R")
source("env_fncs.R")

set.seed(1998)
max_epochs=5000
parameter_mat = make_parameter_mat(max_epochs+100)
set.seed(2001)
total_time=100
num_pulses=4
num_free_pulses=num_pulses-1
state_size = total_time+15
bellmann_error = rep(NA, max_epochs)
waitime_vec=rep(10, num_free_pulses)
state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
potential_actions = 1:15
action_size=length(potential_actions)
actions = rep(NA, max_epochs)
dones = rep(NA, max_epochs)
eps=1
eps.vec=c(eps)
eps_decay=.9995
epsilon_min=.001
minibatch_size=501
epoch=1



rewards = c()
mouse=0
epoch=1

random.init = as.data.frame(matrix(rnorm(length(state_mat)), nrow=nrow(state_mat), ncol=state_size))
random.init$actions=sample(potential_actions,nrow(state_mat), replace=T)
random.init$actions2 = random.init$actions^2

q.fit=readRDS("model.rds")

test_indices=1:100
num_mice = length(test_indices)
optimal_actions=rep(NA, num_mice)
data_mat = matrix(NA, ncol=20, nrow=num_mice)
true.mins = rep(NA, num_mice)
agent.outcomes=rep(NA, num_mice)
day15.outcomes=rep(NA, num_mice)
selected_actions=c()

agent.action.mat = matrix(NA, nrow=num_mice, ncol=num_free_pulses)
reference_days=c(8)

agent.outcome.vec=rep(NA, 100)
refdays = c(10,14,20)
ref.outcome.mat = matrix(NA, nrow=num_mice, ncol=length(refdays)+1)
effect.sizes = rep(NA, ncol(ref.outcome.mat))
references = c(paste("day", refdays), "random")
names(effect.sizes) = references

colnames(ref.outcome.mat) = references

for(mouse in 1:100){
  parameter_vec=parameter_mat[mouse,]
  action.vec=c()
  starting.time=15-2+waitime_vec[1]
  first.action = 15
  one.sequence = generate_one(action.vec, parameter_vec, maxtime=starting.time)
  while(length(action.vec)<num_free_pulses){
  isdone=(length(action.vec)==(num_free_pulses-1))
  one.state=sequence_to_state(one.sequence = as.numeric(one.sequence), action.vec = action.vec, done=isdone, total_time=total_time,num_free_pulses = num_free_pulses)
  one.action = get_action(q.fit,one.state, potential_actions = potential_actions)+waitime_vec[1]
  action.vec=c(action.vec, one.action)
  rt.days = cumsum(action.vec)+15
  seqlength = (max(rt.days)+waitime_vec[1]-2)*(1-isdone)+total_time*isdone
  one.sequence = generate_one(rt.days, parameter_vec, maxtime = seqlength)
  current.time = cumsum(action.vec)+15
  }
  for(ref.idx in 1:length(refdays)){
    ref = refdays[ref.idx]
  one.reference = generate_one(cumsum(rep(ref, num_free_pulses))+15, parameter_vec, maxtime = seqlength)
  ref.outcome.mat[mouse,ref.idx]=log(one.reference[seqlength])
  }
  one.random = sample(potential_actions,num_free_pulses, replace=T)
  one.random = one.random+waitime_vec
  ref.outcome.mat[mouse,length(refdays)+1]= log(generate_one(cumsum(one.random)+15, parameter_vec, maxtime = seqlength))[seqlength]
  agent.action.mat[mouse,]=action.vec
  agent.outcome.vec[mouse] = log(one.sequence[seqlength])
  
}

adj.val=1
plot(density(-agent.outcome.vec+ref.outcome.mat[,1],adjust=adj.val), col="red",type="n", ylim=c(0,0.2),xlab="final ltv reference-final ltv agent", main="4 pulse performance", xlim=c(-.4,.7))
colors= c("red","blue", "orange", "purple")
for(ref.idx in 1:(ncol(ref.outcome.mat))){
  ref = refdays[ref.idx]
  lines(density(ref.outcome.mat[,ref.idx] - agent.outcome.vec, adjust=adj.val+4*(ref.idx==4)), col=colors[ref.idx])
  print(t.test(ref.outcome.mat[,ref.idx] - agent.outcome.vec))
  delta = ref.outcome.mat[,ref.idx] - agent.outcome.vec
  delta2 = ref.outcome.mat[,ref.idx] - ref.outcome.mat[,4]
  effect.sizes[[references[ref.idx]]] = mean(delta)/sd(delta)
  print(summary(delta))
  print(length(delta[delta<0])/100)
  print(sum(delta2<0)/100)
}
legend("topleft", legend=paste(c(paste("day", refdays),"random"), round(effect.sizes, digits=3), sep=": d="), pch=16, col=colors)
abline(v=0, lty=2)
mean(-agent.outcome.vec+ref.outcome.mat[,1])
