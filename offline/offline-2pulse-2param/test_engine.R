
test_indices=1:100
num_mice = length(test_indices)
optimal_actions=rep(NA, num_mice)
data_mat = matrix(NA, ncol=20, nrow=num_mice)
true.mins = rep(NA, num_mice)
agent.outcomes=rep(NA, num_mice)
reference_days=c(8,10,14)
reference.outcomes=matrix(NA, nrow=num_mice, ncol=length(reference_days))
selected_actions=c()

reference_day = 8
for(mouse in 1:num_mice){
  parameter_vec= parameter_mat[mouse,]
  one.reference = generate_one_counterfactualset(parameter_mat[mouse,],total_days = total_days, wait_time = waitime_vec[1], potential_actions = 1:10)
  one.optima= which.min(one.reference[,total_days])+waitime_vec[1]
  optimal_actions[mouse]= as.numeric(one.optima)
  data_mat[mouse,] = one.reference[1,1:20]
  true.mins[mouse]= log(min(one.reference[,total_days]))
  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
  #day15.outcomes[mouse] = one.reference[reference_day-6,total_days]
  reference.outcomes[mouse, refday.idx]=log(one.reference[reference_day-6,total_days])
  }
  one.sequence = data_mat[mouse,]
  one.state=sequence_to_state(one.sequence = as.numeric(one.sequence), action.vec = c(), done=T, num_free_pulses = 1)
  one.action = get_action(q.fit,one.state, potential_actions = 1:10)+7
  selected_actions=c(selected_actions, one.action)
  one.agent.treated.sequence = generate_one(one.action+15, parameter_vec=parameter_vec, maxtime=total_days)
  agent.outcomes[mouse] = log(one.agent.treated.sequence[total_days])
}



selected_actions
deltamodel = abs(optimal_actions-selected_actions)
delta15 =abs(optimal_actions-reference_days[1])
colors = c("red","blue","orange","purple")
plot(density(delta15-deltamodel), xlab="fixed x loss - agent x loss", main="Density plot", col="red",ylim=c(0,0.5), xlim=c(-10,10))
for(refday.idx in 1:length(reference_days)){
  cat("\n info for reference:",reference_days[refday.idx])
  deltaref = abs(optimal_actions-reference_days[refday.idx])
  lines(density(deltaref - deltamodel), col=colors[refday.idx])
  one.delta=deltaref-deltamodel
  print(t.test(deltaref-deltamodel))
  print(100*(length(test_indices)-length(one.delta[one.delta<0]))/length(test_indices))
}
deltarandom = abs(optimal_actions-sample(7:17, num_mice, replace=T))-deltamodel
lines(density(abs(optimal_actions-sample(7:17, num_mice, replace=T)) - deltamodel), col="purple")
legend("topleft", legend=c(paste("day", reference_days)), pch=16, col=colors)
abline(v=0, lty=2)

deltamodel = agent.outcomes-true.mins
delta15 =reference.outcomes[,1]-true.mins
colors = c("red","blue","orange","purple")
plot(density(delta15-deltamodel), xlab="fixed y loss - agent y loss", main="Density plot", col="red", ylim=c(0,60))
for(refday.idx in 2:length(reference_days)){
  deltaref = reference.outcomes[,refday.idx]-true.mins
  lines(density(deltaref - deltamodel), col=colors[refday.idx])
  cat(deltamodel[deltamodel>deltaref])
}
lines(density(abs(optimal_actions-sample(7:17, num_mice, replace=T)) - deltamodel), col="purple")
legend("topleft", legend=c(paste("day", reference_days), "random"), pch=16, col=c(colors, "purple"))
abline(v=0, lty=2)

table(selected_actions, optimal_actions)
