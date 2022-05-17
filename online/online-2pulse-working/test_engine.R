setwd("C:/Users/s198663/Documents/RPAI/online/online-2pulse-working")
source("data_generation_fncs.R")
source("env_fncs_2pulse.R")
set.seed(1998)
test_num=50
maxtime = 40
num_free_pulses=1
wait_time=9
potential_actions = 1:5 + wait_time
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
inc_days = 15-2+wait_time
minibatch_parameters=make_parameter_mat(test_num)

agent.outcomes=rep(NA, test_num)
reference_days=c(1,potential_actions,20)
reference.outcomes=matrix(NA, nrow=test_num, ncol=length(reference_days)+1)
colnames(reference.outcomes) = c(paste("day",reference_days),'random')
selected_actions=matrix(NA, nrow=test_num, ncol=num_free_pulses)
est.par.mat = matrix(NA, nrow=test_num, ncol=3)
colnames(est.par.mat) = c("mu","lambda","rho")
optima.verif = rep(NA, test_num)
optimal_actions = rep(NA,test_num)

for(mouse in 1:test_num){
  parameter_vec= minibatch_parameters[mouse,]
  one.reference = generate_one_counterfactualset(minibatch_parameters[mouse,],total_days = maxtime, wait_time = 0, potential_actions = 1:20)
  one.optima= which.min(one.reference[,maxtime])+wait_time
  optimal_actions[mouse] = one.optima
  optima.verif[mouse] = get.optim.plan(parameter_vec=parameter_vec, maxtime=maxtime, all.action.mat = all.action.mat)
  optimal_actions[mouse]= as.numeric(one.optima)

  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
    #day15.outcomes[mouse] = one.reference[reference_day-6,total_days]
    reference.outcomes[mouse, refday.idx]=log(one.reference[reference_day,maxtime])
  }
  
  random.action = sample(potential_actions, size=num_free_pulses, replace=T)
  reference.outcomes[mouse,'random'] = log(one.reference[random.action-wait_time,maxtime])
  one.sequence = generate_one(radiation_days = c(), parameter_vec = parameter_vec, maxtime=inc_days)$ltv
  estimated.par = optimize.model(one.sequence)
  names(estimated.par) = c("mu","lambda","rho")
  est.par.mat[mouse,] = estimated.par
  if(mouse>1){weight.factor=apply(est.par.mat,c(2),sd, na.rm=T)*c(150,3,5)}else{weight.factor=1}
  predictive.params = estimated.par*weight.factor+colMeans(est.par.mat, na.rm=T)*(1-weight.factor)
  #print(estimated.par)
  one.action = get.optim.plan(parameter_vec=predictive.params,maxtime=maxtime, all.action.mat = all.action.mat)
  selected_actions[mouse,] = one.action
  one.agent.treated.sequence = generate_one(cumsum(one.action)+15, parameter_vec=parameter_vec, maxtime=maxtime+10)$ltv
  agent.outcomes[mouse] = log(one.agent.treated.sequence[maxtime])
}


reference_actions = c(reference_days, "random")

effect.sizes = rep(NA, 4)

#x11()
results = matrix(NA, nrow=ncol(reference.outcomes),ncol=4)
colnames(results) = c("means","medians", "mins", "maxes")
rownames(results) = colnames(reference.outcomes)
par(mfrow=c(2,2))
colors = c("red","blue","orange","purple")
plot(density(agent.outcomes-reference.outcomes[,1]),xlim=c(-.05,.2),xlab="fixed final ltv - agent final ltv",type="n", main="2 pulse performance", col="red", ylim=c(0,60))
for(refday.idx in 1:ncol(reference.outcomes)){
  deltaref = reference.outcomes[,refday.idx]
  lines(density(deltaref - agent.outcomes), col=colors[refday.idx])
  test=t.test(deltaref-agent.outcomes)
  print(test)
  loss=deltaref-agent.outcomes
  cohen = mean(deltaref-agent.outcomes)/sd(deltaref-agent.outcomes)
  cat("effect size: ", cohen)
  effect.sizes[refday.idx] = cohen
  results[refday.idx,] = c(mean(loss), median(loss), min(loss), max(loss))
  #cat(agent.outcomes[agent.outcomes>deltaref])
}

#effect.sizes[4] = mean(delta.diff.random)/sd(delta.diff.random)
#lines(density(abs(optimal_actions-sample(potential_actions, test_num, replace=T)) - agent.outcomes), col="purple")
legend("topright", legend=paste(c(paste("day", reference_days), "random"), round(effect.sizes, 3), sep=": d="), pch=16, col=c(colors, "purple"))
abline(v=0, lty=2)

par(mfrow=c(1,3))
plot(minibatch_parameters[,'lambda'], est.par.mat[,'lambda'],pch=16, main="lambda", ylab="estimated", xlab="actual")
plot(minibatch_parameters[,'rho'], est.par.mat[,'rho'], pch=16,main="rho", ylab="estimated", xlab="actual")
plot(minibatch_parameters[,'mu'], est.par.mat[,'mu'], pch=16,main="mu", ylab="estimated", xlab="actual")
