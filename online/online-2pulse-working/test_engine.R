test_indices=1:test_num

optimal_actions=rep(NA, test_num)
data_mat = matrix(NA, ncol=20, nrow=test_num)
true.mins = rep(NA, test_num)
agent.outcomes=rep(NA, test_num)
reference_days=c(10,14,20)
reference.outcomes=matrix(NA, nrow=test_num, ncol=length(reference_days))
selected_actions=c()

reference_day = 8
for(mouse in 1:test_num){
  parameter_vec= parameter_mat[mouse,]
  one.reference = generate_one_counterfactualset(parameter_mat[mouse,],total_days = total_days+5, wait_time = waitime_vec[1], potential_actions = 1:15)
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
  one.action = get_action(q.fit,one.state, potential_actions = 1:13)+waitime_vec[1]
  selected_actions=c(selected_actions, one.action)
  one.agent.treated.sequence = generate_one(one.action+15, parameter_vec=parameter_vec, maxtime=total_days)
  agent.outcomes[mouse] = log(one.agent.treated.sequence[total_days])
}


reference_actions = c(reference_days, "random")
selected_actions
deltamodel = abs(optimal_actions-selected_actions)
delta15 =abs(optimal_actions-reference_days[1])
colors = c("red","blue","orange","purple")
adj_factor=2
plot(density(delta15-deltamodel, adjust = adj_factor), xlab="fixed x loss - agent x loss", main="Density plot", col="red",ylim=c(0,0.8), xlim=c(-10,10))
for(refday.idx in 1:length(reference_days)){
  cat("\n info for reference:",reference_days[refday.idx])
  deltaref = abs(optimal_actions-reference_days[refday.idx])
  lines(density(deltaref - deltamodel, adjust = adj_factor), col=colors[refday.idx])
  one.delta=deltaref-deltamodel
  print(t.test(deltaref-deltamodel))
  cohen = mean(deltaref-deltamodel)/sd(deltaref-deltamodel)
  cat("effect size: ", cohen)
  print(100*(length(test_indices)-length(one.delta[one.delta<0]))/length(test_indices))
}
deltarandom = abs(optimal_actions-sample(7:17, test_num, replace=T))-deltamodel
lines(density(abs(optimal_actions-sample(7:17, test_num, replace=T)) - deltamodel, adjust=adj_factor), col="purple")
legend("topleft", legend=c(paste("day", reference_days),"random"), pch=16, col=colors)
abline(v=0, lty=2)
effect.sizes = rep(NA, 4)

deltamodel = agent.outcomes-true.mins
delta15 =reference.outcomes[,1]-true.mins
colors = c("red","blue","orange","purple")
plot(density(delta15-deltamodel), xlab="fixed personalization loss - agent personalization loss",type="n",xlim=c(-.2,.2), main="2 pulse performance", col="red", ylim=c(0,60))
for(refday.idx in 1:3){
  deltaref = reference.outcomes[,refday.idx]-true.mins
  lines(density(deltaref - deltamodel), col=colors[refday.idx])
  test=t.test(deltaref-deltamodel)
  print(test)
  cohen = mean(deltaref-deltamodel)/sd(deltaref-deltamodel)
  cat("effect size: ", cohen)
  effect.sizes[refday.idx] = cohen
  #cat(deltamodel[deltamodel>deltaref])
}
delta.diff.random = abs(optimal_actions-sample(7:17, test_num, replace=T)) - deltamodel
effect.sizes[4] = mean(delta.diff.random)/sd(delta.diff.random)
lines(density(abs(optimal_actions-sample(7:17, test_num, replace=T)) - deltamodel), col="purple")
legend("topright", legend=paste(c(paste("day", reference_days), "random"), round(effect.sizes, 3), sep=": d="), pch=16, col=c(colors, "purple"))
abline(v=0, lty=2)

table(selected_actions, optimal_actions)
