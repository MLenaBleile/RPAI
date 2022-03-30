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

agent.outcome.vec=rep(NA, num_mice)
refdays = c(10,14,20)
ref.outcome.mat = matrix(NA, nrow=num_mice, ncol=length(refdays)+1)
effect.sizes = rep(NA, ncol(ref.outcome.mat))
references = c(paste("day", refdays), "random")
names(effect.sizes) = references

random.action.mat = agent.action.mat
colnames(ref.outcome.mat) = references

for(mouse in 1:num_mice){
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
    rt.days = cumsum(action.vec)+first.action
    seqlength = (max(rt.days)+waitime_vec[1]-2)*(1-isdone)+total_time*isdone
    one.sequence = generate_one(rt.days, parameter_vec, maxtime = seqlength)
    current.time = cumsum(action.vec)+first.action
  }
  for(ref.idx in 1:length(refdays)){
    ref = refdays[ref.idx]
    one.reference = generate_one(cumsum(c(ref,ref, ref))+first.action, parameter_vec, maxtime = seqlength)
    ref.outcome.mat[mouse,ref.idx]=log(one.reference[seqlength])
  }
  one.random = sample(potential_actions,num_free_pulses, replace=T)
  random.action.mat[mouse,] = one.random
  one.random = one.random+waitime_vec
  ref.outcome.mat[mouse,length(refdays)+1]= log(generate_one(cumsum(one.random)+15, parameter_vec, maxtime = seqlength))[seqlength]
  agent.action.mat[mouse,]=action.vec
  agent.outcome.vec[mouse] = log(one.sequence[seqlength])
  
}

adj.val=1
plot(density(-agent.outcome.vec+ref.outcome.mat[,1],adjust=adj.val), col="red",type="n", ylim=c(0,5),xlab="final ltv reference-final ltv agent", main="4 pulse performance", xlim=c(-.4,.7))
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
legend("topright", legend=paste(c(paste("day", refdays),"random"), round(effect.sizes, digits=3), sep=": d="), pch=16, col=colors)
abline(v=0, lty=2)
mean(-agent.outcome.vec+ref.outcome.mat[,1])
