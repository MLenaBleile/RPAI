setwd("C:/Users/s198663/Documents/RPAI/online/online-2pulse-working")
source("datagen_fncs.R")
source("env_fncs_2pulse.R")
library(brms)
library(nlme)
library(boot)
set.seed(1998)
calibration_num=100
test_num=calibration_num+190
recompute_marginals=T
mod.to.fit="sumexp"
gen.mod="sumexp"
loss.weights=c(1,0,0)
names(loss.weights) = c("raw","liklihood","effects")
###parameters that we don't assume to be known
param_names=c("rho1","theta1","beta0","beta1","gg")
random_params=c('gg')
fixed_params = setdiff(param_names, random_params)
optimal_actions = rep(NA,test_num)


maxtime = 40
num_free_pulses=1
wait_time=7
potential_actions = 1:7 + wait_time
all.action.mat = make_potential_action_mat(potential_actions,1)
bounds = matrix(NA,nrow=length(param_names), ncol=2)
colnames(bounds) = c("upper","lower")
rownames(bounds) = param_names
bounds[,'upper'] = 1
bounds[,'lower'] = 0
bounds['theta1','upper']=.3


inc_days = 15+wait_time
minibatch_parameters=make_parameter_mat(test_num, gen.mod)

agent.outcomes=rep(NA, test_num)
reference_days=c(1,potential_actions,20)
reference.outcomes=matrix(NA, nrow=test_num, ncol=length(reference_days)+1)
colnames(reference.outcomes) = c(paste("day",reference_days),'random')
selected_actions=matrix(NA, nrow=test_num, ncol=num_free_pulses)
all.predictive.pars = matrix(NA, nrow=test_num, ncol=length(param_names))
colnames(all.predictive.pars) = param_names

all.pairs = array(dim=c(test_num,maxtime,4))
dimnames(all.pairs) = list(NULL, NULL, c("ltv","d","p",'day'))
names(dimnames(all.pairs)) = c('animal','time', 'feature')



agent.actions=numeric(test_num)
mouse=1
estimated.par = NULL

for(mouse in mouse:test_num){
  cat("optimizing mouse ",mouse,"\n")
  true.param.vec= minibatch_parameters[mouse,]
  one.reference.pairs = generate_one_counterfactualset(parameter_vec = true.param.vec,total_days = maxtime,potential_actions = 1:20, gen.mod=gen.mod)
  ##first check where the true optimum is and calculate each reference curve
  one.optima= which.min(one.reference.pairs[potential_actions,maxtime,'ltv'])
  optimal_actions[mouse] = potential_actions[one.optima]
  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
    reference.outcomes[mouse, refday.idx]=one.reference.pairs[reference_day,maxtime,'ltv']
  }
  
  random.action = sample(potential_actions, size=num_free_pulses, replace=T)
  reference.outcomes[mouse,'random'] = one.reference.pairs[random.action,maxtime,'ltv']
  inc.pair = generate_one(radiation_days=c(), parameter_vec=true.param.vec, maxtime=inc_days,add.measurement.noise=T, gen.mod=gen.mod)
  all.pairs[mouse,1:inc_days,] = inc.pair[1,1:inc_days,]
  ###if in the calibration phase, generate a random action
  if(mouse<calibration_num){
    one.action = sample(potential_actions,size=num_free_pulses, replace=T)
    agent.actions[mouse]=-9
  }else{
    
    ##melt the all.pairs data 
    all.ltv = all.pairs[1:mouse,,'ltv']
    dimnames(all.ltv) = list(c(1:mouse), 1:maxtime)
    names(dimnames(all.ltv)) = c("animal","current.time")
    all.pairs.molten = reshape2::melt(all.ltv, value.name="ltv")
    all.pairs.molten$d1 = apply(all.pairs[1:mouse,,'d'],c(1), get.chg.idx)
    #all.pairs.molten$ltv[all.pairs.molten$animal==mouse& all.pairs.molten$current.time<= inc_days]=inc.pair[1,]
    all.pairs.clean=all.pairs.molten[all.pairs.molten$current.time>1&(!is.na(all.pairs.molten$ltv)),]
    
    estimated.par = optimize.model(all.pairs[1:(mouse),,], param_names, fixed_param_vec = c(), loss.weights=loss.weights,mod.to.fit="sumexp", bounds=bounds, estimated.par=estimated.par)
    get.mean.sumexp = Vectorize(function(beta0,beta1,gg, d1,current.time){
      scale.fact=1
      d0=sqrt(15)/scale.fact
      d1 = sqrt(d1)/scale.fact
      rho1 = estimated.par['rho1']
      rhoo = estimated.par['rho1']
      beta0 = estimated.par['beta0']

      theta1 = estimated.par['theta1']
      rho2=rho1*(1- exp(-abs(theta1)*(d1-d0)))
      gg=abs(gg)
      #print(length(d1))
      #print(rho2)
      time.vec=sqrt(1:maxtime)/scale.fact
      first.decay.term = rhoo*exp(beta0+beta1*d0)*exp(-gg*(time.vec-d0))*(time.vec>d0)
      fct.mean.exp = first.decay.term + (1-rhoo*(time.vec>d0))*exp(beta0+beta1*time.vec)
      if(maxtime>d1){
        fct.mean.exp.1=fct.mean.exp
        fct.mean.exp = rho2*(time.vec>d1)*fct.mean.exp.1*exp(-gg*(time.vec-d1))+(1-rho2*(time.vec>d1))*fct.mean.exp.1
        
      }else{
        
      }
      out=log(fct.mean.exp)[maxtime]
      return(as.numeric(out))
    }, vectorize.args = c("beta0","beta1","gg", "d1","current.time"))
    fixed_param_vec = estimated.par[fixed_params]
    one.mouse = array(all.pairs[mouse,!(is.na(all.pairs[mouse,,'ltv'])),], dim = c(1, dim(all.pairs[mouse,,])), dimnames=dimnames(all.pairs))

    indiv.estimated.par = optimize.model(one.mouse, param_names,fixed_param_vec=fixed_param_vec, loss.weights=loss.weights, bounds=bounds,mod.to.fit="sumexp")
    predictive.par=estimated.par
    predictive.par[random_params] = indiv.estimated.par[random_params]
    all.predictive.pars[mouse,param_names] = predictive.par[param_names] 
    one.action = get.optim.plan(predictive.par,maxtime, all.action.mat, "sumexp")
    agent.actions[mouse]=one.action
  }
  
  one.agent.treated.sequence = generate_one(cumsum(one.action)+15, parameter_vec=true.param.vec, maxtime=maxtime, gen.mod=gen.mod, add.measurement.noise = T)
  all.pairs[mouse,-c(1:inc_days),] = one.agent.treated.sequence[1,-c(1:inc_days),]
  one.noiseless.sequence = generate_one(cumsum(one.action)+15, parameter_vec=true.param.vec, maxtime=maxtime, gen.mod=gen.mod, add.measurement.noise = F)
  agent.outcomes[mouse] = one.noiseless.sequence[1,maxtime,'ltv']
}


###check how the agent did as compared to each reference
results = matrix(NA, nrow=ncol(reference.outcomes),ncol=5)
colnames(results) = c("mean perc reduc","median perc reduc", "midpoint", "prop better","prop asgoodas")
rownames(results) = colnames(reference.outcomes)

for(refday.idx in 1:ncol(reference.outcomes)){
  deltaref = reference.outcomes[,refday.idx]
  loss=deltaref-agent.outcomes
  loss=loss[71:103]
  plot(density(loss), main=colnames(reference.outcomes)[refday.idx])
  results[refday.idx,] = c(mean(loss), median(loss), (max(loss)-abs(min(loss))), sum(loss>0)/length(loss[loss!=0]),sum(loss>=0)/length(loss) )
}
print(results)


beepr::beep(sound=0)


