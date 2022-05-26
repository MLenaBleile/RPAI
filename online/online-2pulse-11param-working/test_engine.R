setwd("C:/Users/s198663/Documents/RPAI/online/online-2pulse-11param-working")
source("data_generation_fncs.R")
source("env_fncs_2pulse.R")
set.seed(1998)
calibration_num=10
test_num=calibration_num+45
recompute_marginals=T
changepoint.time= calibration_num+1
maxtime = 40
num_free_pulses=1
wait_time=9
potential_actions = 1:5 + wait_time
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
##this is the date of the first IO treatment (two days before the first RT), plus the waiting time.
inc_days = 15-2+wait_time
minibatch_parameters=make_parameter_mat(test_num)


##initialize a bunch of objects for storing the outcomes
agent.outcomes=rep(NA, test_num)
reference_days=c(1,potential_actions,20)
reference.outcomes=matrix(NA, nrow=test_num, ncol=length(reference_days)+1)
colnames(reference.outcomes) = c(paste("day",reference_days),'random')
selected_actions=matrix(NA, nrow=test_num, ncol=num_free_pulses)
###parameters that we don't assume to be known
param_names = c("mu","lambda","rho","phi","alpha1","gam",'tau','omega1','omega2','alpha','beta')
est.par.mat = matrix(NA, nrow=test_num, ncol=length(param_names))
colnames(est.par.mat) = param_names
predictive.param.mat = est.par.mat
optima.verif = rep(NA, test_num)
optimal_actions = rep(NA,test_num)
##parameters that we assume have a distribution are the random_parameters
random_params=c('mu','lambda','rho')
fixed_params = setdiff(param_names, random_params)
weight.factors = matrix(nrow=test_num, ncol=length(random_params))
initial.sequences = matrix(nrow=test_num, ncol=inc_days)
all.pairs = array(dim=c(test_num,maxtime,4))
dimnames(all.pairs) = list(NULL, NULL, c("ltv","d","p",'day'))
names(all.pairs) = c('animal','time', 'feature')

if(recompute_marginals){
  overall.estimate.mat = matrix(NA, nrow=test_num, ncol=length(param_names)+1)
  colnames(overall.estimate.mat) = c(param_names,"sigma")
}else{}

all.counterfactual.pairs = array(dim=c(0,maxtime,4))
start.idx = 1
#if(mouse<test_num){start.idx=mouse}

for(mouse in start.idx:test_num){
  cat("optimizing mouse ",mouse,"\n")
  parameter_vec= minibatch_parameters[mouse,]
  one.reference.pairs = generate_one_counterfactualset(minibatch_parameters[mouse,],total_days = maxtime,potential_actions = 1:20)
  one.optima= which.min(one.reference.pairs[potential_actions,maxtime,'ltv'])
  optimal_actions[mouse] = one.optima+wait_time
  optima.verif[mouse] = get.optim.plan(parameter_vec=parameter_vec, maxtime=maxtime, all.action.mat = all.action.mat)
  

  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
    #day15.outcomes[mouse] = one.reference[reference_day-6,total_days]
    reference.outcomes[mouse, refday.idx]=one.reference.pairs[reference_day,maxtime,'ltv']
  }
  
  random.action = sample(potential_actions, size=num_free_pulses, replace=T)
  reference.outcomes[mouse,'random'] = one.reference.pairs[random.action,maxtime,'ltv']
  inc.pair = generate_one(radiation_days = c(), parameter_vec = parameter_vec, maxtime=inc_days)
  initial.sequences[mouse,] = inc.pair[1,,'ltv']
  estimated.par = optimize.model(pair.set=inc.pair, param_names = param_names, fixed_param_vec=c())
  names(estimated.par) = param_names
  est.par.mat[mouse,fixed_params] = estimated.par[fixed_params]
  predictive.params = estimated.par
  if(mouse>2){
    mssdiff = sqrt(sum((inc.pair[1,,'ltv']-colMeans(initial.sequences, na.rm=T))^2))
    weight.factor=1+sqrt(apply(est.par.mat,c(2),sd, na.rm=T)*mssdiff)
    if(recompute_marginals){
      overall.params= optimize.model(pair.set=all.pairs[1:(mouse-1),,], param_names = param_names, fixed_param_vec=c(), liklihood.based = (mouse>=changepoint.time),starting.params=colMeans(overall.estimate.mat, na.rm=T))
      overall.estimate.mat[mouse,names(overall.params)]= overall.params
      write.csv(overall.estimate.mat, "overall.estimate.mat.csv")
    }else{overall.params = overall.estimate.mat[mouse,]}
    
  fixed.weight=.5*apply(overall.estimate.mat[,fixed_params],c(2), sd, na.rm=T)/apply(est.par.mat[,fixed_params],c(2), sd, na.rm=T)
  fixed.weight[fixed.weight>1]=1
  predictive.params[fixed_params] = apply(est.par.mat[,fixed_params],c(2), median, na.rm=T)*fixed.weight +colMeans(overall.estimate.mat, na.rm=T)[fixed_params]*(1-fixed.weight)
  
  }else{
    weight.factor=rep(1, length(param_names))
    names(weight.factor)=param_names
    overall.estimate.mat[mouse,names(estimated.par)] = estimated.par
  }
  weight.factor[weight.factor>1]=1
  #weight.factor=rep(1, length(param_names))
  #names(weight.factor) = param_names
  weight.factors[mouse,] = weight.factor[random_params]
  weight.factor.random=weight.factor[random_params]
  fixed_param_vec=colMeans(overall.estimate.mat[,fixed_params], na.rm=T)
  reestimated.par = optimize.model(pair.set=inc.pair, param_names = random_params,fixed_param_vec=fixed_param_vec)
  est.par.mat[mouse,random_params] = reestimated.par[random_params]
  predictive.params[random_params] = reestimated.par[random_params]*weight.factor.random+colMeans(est.par.mat[,random_params], na.rm=T)*(1-weight.factor.random)
  
 
  predictive.param.mat[mouse,]=predictive.params
 
  if(mouse<calibration_num){
    one.action=sample(potential_actions, size=num_free_pulses, replace=T)
  }else{
    one.action = get.optim.plan(parameter_vec=predictive.params,maxtime=maxtime, all.action.mat = all.action.mat)
    selected_actions[mouse,] = one.action
  }
  selected_actions[mouse,] = one.action
  one.agent.treated.sequence = generate_one(cumsum(one.action)+15, parameter_vec=parameter_vec, maxtime=maxtime)
  all.pairs[mouse,,] = one.agent.treated.sequence[1,,]
  one.cf.reference.pairs = generate_one_counterfactualset(predictive.params,total_days = maxtime,potential_actions = potential_actions)
  one.cf.reference.pairs[one.action-wait_time,,] = one.agent.treated.sequence[1,,]
  all.counterfactual.pairs= abind::abind(all.counterfactual.pairs, one.cf.reference.pairs, along=1)
  agent.outcomes[mouse] = one.agent.treated.sequence[1,maxtime,'ltv']
}


reference_actions = c(reference_days, "random")

effect.sizes = rep(NA, 5)

#x11()
results = matrix(NA, nrow=ncol(reference.outcomes),ncol=5)
colnames(results) = c("mean","median", "midpoint", "prop better","prop asgoodas")
rownames(results) = colnames(reference.outcomes)
par(mfrow=c(2,2))
colors = c("red","blue","orange","purple")
plot(density(agent.outcomes-reference.outcomes[,1]),xlim=c(-.05,.2),xlab="fixed final ltv - agent final ltv",type="n", main="2 pulse performance", col="red", ylim=c(0,60))
for(refday.idx in 1:ncol(reference.outcomes)){
  deltaref = reference.outcomes[,refday.idx]
  #lines(density(deltaref[!is.na(deltaref)] - agent.outcomes[!is.na(deltaref)]), col=colors[refday.idx])
  loss=deltaref-agent.outcomes
  loss=loss[-c(1:calibration_num)]
  results[refday.idx,] = c(mean(loss, na.rm=T), median(loss, na.rm=T), max(loss)-min(loss), sum(loss>0)/sum(loss!=0),sum(loss>=0)/length(loss) )
  #cat(agent.outcomes[agent.outcomes>deltaref])
}

beepr::beep(sound=8)

#effect.sizes[4] = mean(delta.diff.random)/sd(delta.diff.random)
#lines(density(abs(optimal_actions-sample(potential_actions, test_num, replace=T)) - agent.outcomes), col="purple")
legend("topright", legend=paste(c(paste("day", reference_days), "random"), round(effect.sizes, 3), sep=": d="), pch=16, col=c(colors, "purple"))
abline(v=0, lty=2)

par(mfrow=c(2,2))
#plot(minibatch_parameters[,'lambda'], est.par.mat[,'lambda'],pch=16, main="lambda", ylab="estimated", xlab="actual")
#plot(minibatch_parameters[,'rho'], est.par.mat[,'rho'], pch=16,main="rho", ylab="estimated", xlab="actual")
#plot(minibatch_parameters[,'mu'], est.par.mat[,'mu'], pch=16,main="mu", ylab="estimated", xlab="actual")
colnames(predictive.param.mat)=colnames(est.par.mat)
colnames(weight.factors) = random_params
for(param_name in random_params){
  plot(minibatch_parameters[,param_name], est.par.mat[, param_name], pch=16,ylab="direct estimate", main=param_name)
  plot(minibatch_parameters[,param_name], predictive.param.mat[, param_name], pch=16,ylab="predictive", main=param_name)
  plot(minibatch_parameters[,param_name], overall.estimate.mat[, param_name], pch=16,ylab="overall", main=param_name)
  hist(weight.factors[,param_name])
  }
print(results)

fixed.weights=matrix(nrow=test_num, ncol=length(fixed_params))
colnames(fixed.weights)=fixed_params
for(mm in 2:test_num){
  fixed.weights[mm, fixed_params] = apply(overall.estimate.mat[1:mm,fixed_params],c(2), sd, na.rm=T)/apply(est.par.mat[1:mm,fixed_params],c(2), sd, na.rm=T)
}