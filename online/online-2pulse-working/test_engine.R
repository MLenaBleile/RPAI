test_indices=1:50
num_mice = length(test_indices)
optimal_actions=rep(NA, num_mice)
data_mat = matrix(NA, ncol=20, nrow=num_mice)
true.mins = rep(NA, num_mice)
agent.outcomes=rep(NA, num_mice)
reference_days=c(1,10:14,20)
reference.outcomes=matrix(NA, nrow=num_mice, ncol=length(reference_days)+1)
references=c(paste("day", reference_days),"random")
colnames(reference.outcomes) = references
selected_actions=c()


for(mouse in 1:num_mice){
  parameter_vec= parameter_mat[mouse,]
  one.reference = generate_one_counterfactualset(parameter_mat[mouse,],total_days = total_days+5, wait_time = 0, potential_actions = 1:20)
  one.optima= which.min(one.reference[,total_days])
  optimal_actions[mouse]= as.numeric(one.optima)
  #data_mat[mouse,] = one.reference[1,1:20]
  true.mins[mouse]= log(min(one.reference[,total_days]))
  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
    #day15.outcomes[mouse] = one.reference[reference_day-6,total_days]
    reference.outcomes[mouse, refday.idx]=log(one.reference[reference_day,total_days])
  }
  one.random.action=sample(potential_actions+waitime_vec[1],1)
  reference.outcomes[mouse,'random'] = log(one.reference[one.random.action,total_days])
  current.time=15-2+waitime_vec[1]
  one.sequence = generate_one(c(), parameter_vec, current.time)
  one.state=sequence_to_state(one.sequence = as.numeric(one.sequence), action.vec = c(), done=T, num_free_pulses = 1)
  one.action = get_action(q.fit,one.state, potential_actions = potential_actions)+waitime_vec[1]
  selected_actions=c(selected_actions, one.action)
  one.agent.treated.sequence = generate_one(one.action+15, parameter_vec=parameter_vec, maxtime=total_days)
  agent.outcomes[mouse] = log(one.reference[one.action,total_days])
}

overall.min = min(c(agent.outcomes, as.numeric(reference.outcomes)))
overall.max = max(c(agent.outcomes, as.numeric(reference.outcomes)))
result = matrix(NA, nrow=ncol(reference.outcomes), ncol=4)
colnames(result) = c("mean","median","max", "min")
rownames(result) = colnames(reference.outcomes)
perc.better=c()
perc.asgoodas =c()
for(jj in 1:length(references)){
  loss = agent.outcomes- reference.outcomes[,jj] 
  result[jj,] = c(mean(exp(loss)), median(exp(loss)), min(exp(loss)), max(exp(loss)))
  perc.better[jj] = sum(loss<0)/sum(loss!=0)
  perc.asgoodas[jj] = sum(loss<=0)/length(loss)
  }

print(result)


# deltamodel = agent.outcomes - reference.outcomes[,'day 1']
# 
# colors = c("red","blue","orange","purple")
# plot(density(delta15-deltamodel), xlab="fixed personalization loss - agent personalization loss",type="n",xlim=c(-.2,.5), main="2 pulse performance", col="red", ylim=c(0,60))
# for(refday.idx in 1:ncol(reference.outcomes)){
#   deltaref = reference.outcomes[,refday.idx]-reference.outcomes[,'day 1']
#   lines(density(exp(deltaref - deltamodel)), col=colors[refday.idx])
#   test=t.test(deltaref-deltamodel)
#   print(test)
#   cohen = mean(deltaref-deltamodel)
#   cat("effect size: ", cohen)
#   cat("minimum:",min(deltaref-deltamodel))
#   effect.sizes[refday.idx] = cohen
#   #cat(deltamodel[deltamodel>deltaref])
# }
# #delta.diff.random = abs(optimal_actions-sample(7:17, num_mice, replace=T)) - deltamodel
# names(effect.sizes) = colnames(reference.outcomes)
# #lines(density(abs(optimal_actions-sample(7:17, num_mice, replace=T)) - deltamodel), col="purple")
# legend("topright", legend=paste(c(paste("day", reference_days), "random"), round(effect.sizes, 3), sep=": "), pch=16, col=c(colors, "purple"))
# abline(v=0, lty=2)

table(selected_actions, optimal_actions)
