test_indices=1:test_num
optimal_actions=rep(NA, test_num)
true.mins = rep(NA, test_num)
agent.outcomes=rep(NA, test_num)
day15.outcomes=rep(NA, test_num)
selected_actions=c()

agent.action.mat = matrix(NA, nrow=test_num, ncol=num_free_pulses)

agent.outcome.vec=rep(NA, test_num)
refdays = c(1, potential_actions+wait_time,20)
ref.outcome.mat = matrix(NA, nrow=test_num, ncol=length(refdays)+1)
references=c(paste("day", refdays), "random")
colnames(ref.outcome.mat) = c(paste("day", refdays), "random")

for(mouse in 1:test_num){
  parameter_vec=parameter_mat[mouse+max_epochs/num_free_pulses,]
  action.vec=c()
  one.sequence = generate_one(action.vec, parameter_vec, maxtime=inc.days)
  while(length(action.vec)<num_free_pulses){
    isdone=(length(action.vec)==(num_free_pulses-1))
    one.state=sequence_to_state(one.sequence = as.numeric(one.sequence), action.vec = action.vec, wait_time=wait_time,total_time=total_time,done=isdone, num_free_pulses = 2)
    one.action = get_action(q.fit,one.state, potential_actions = potential_actions)+wait_time
    action.vec=c(action.vec, one.action)
    seqlength = (inc.days+sum(action.vec))*(1-isdone)+total_time*isdone
    one.sequence = generate_one(cumsum(action.vec)+15, parameter_vec, maxtime = seqlength)
  }
  for(ref.idx in 1:length(refdays)){
    ref = refdays[ref.idx]
    one.reference = generate_one(cumsum(c(ref,ref))+15, parameter_vec, maxtime = seqlength)
    ref.outcome.mat[mouse,ref.idx]=log(one.reference[seqlength])
  }
  one.random = sample(potential_actions+wait_time,num_free_pulses, replace=T)
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
  print(1-length(delta[delta<=0])/test_num)
  cat("proportion of times agent at least as good as", colnames(ref.outcome.mat)[ref.idx],"day spacing:" )
  print(1-length(delta[delta<0])/test_num)
  cat("proportion of times ", colnames(ref.outcome.mat)[ref.idx],"was better than random" )
  print(1-length(delta[delta2<=0])/test_num)
}
legend("topright", legend=paste(c(paste("day", refdays),"random"), round(effect.sizes, digits=3), sep=": d="), pch=16, col=colors)
abline(v=0, lty=2)
print(exp(mean(agent.outcome.vec-ref.outcome.mat[,1])))
