make_potential_action_mat = function(potential_actions, num_free_pulses){
  expand.grid(rep(list(potential_actions), num_free_pulses))
}

get.optim.plan = function(fit1, mouse, potential_actions, maxtime){
  predicted.outcomes=numeric(length(potential_actions))
  names(predicted.outcomes) = as.character(potential_actions)
  for(action in potential_actions){
    one.cf.input = matrix(c(mouse, maxtime, action), nrow = 1)
    colnames(one.cf.input) = c("animal", "time", "d1")
    predicted.outcomes[as.character(action)]=posterior_predict(fit1, newdata=as.data.frame(one.cf.input))
    
  }
  potential_actions[which.min(predicted.outcomes)]
}

get.chg.idx = function(x){
  chgpts=which(x[-1] != x[-length(x)])+1
  d1.pt=chgpts[chgpts!=15]
  if(length(d1.pt)==0){
    ##assign random action (or Infty) to the incomplete sequence
    d1.pt=80
  }
  d1.pt
}

make_default_para_set = function(){
  para_set = matrix(NA, nrow=27, ncol=3)
  param_names = c("alpha","beta","omega1","omega2","phi","lambda","rho","BT0","alpha1","TT0","psi","mu","gam","tau")
  group_names = c("Vehicle","PDL1","10GyD01","10GyD01+PDL1","10GyD010","10GyD010+PDL1","10GyD011011","10GyD011011+PDL1","10GyD020","10GyD020+PDL1","10GyD012021","10GyD012021+PDL1","Tinit")
  rownames(para_set) = c(param_names, group_names)
  colnames(para_set) = c("lower","upper","best")
  para_set[,"lower"] = c(rep(0,13),rep(1,14))
  para_set[,"upper"] = 1
  para_set[c(7,8,10),"upper"]=3
  para_set[15:27,"upper"] = 10
  best_penelope_params = c(0.02404,0.00148,0.02349,0.01404,0.96396,0.30441,1.70781,0.00001,0.05205,0.00000,0.37045,0.21601,0.88301,8.79821,10.0000,9.99999,1.27574,2.12489,3.12656,1.21461,1.72499,3.49761,4.53273,4.23883,3.18113,6.42940,1)
  default_best_params = c(NA,NA,.0029,0.001,1, 0.387,0.441,6e-6,0.015,.000,NA, 0.231,NA,NA,1.452,8.784)
  default_best_group_params = c(9.460,10,1.643,2.100,3.665,0.657,2.333,2.568,3.143,3.288,6.615,2.685)
  para_set[c(param_names, group_names),'best'] = best_penelope_params
  #para_set[c("D0","dq"),"best"] = c(1.452,8.784)
  para_set
}
truncnorm = Vectorize(function(samples, loc=0, scale=1, upr=Inf, lwr=-Inf){
  proposed = rnorm(samples*2, mean=loc,sd=sqrt(scale))
  if(is.nan(proposed[1])){
    cat("scale:",sqrt(scale), "loc:",loc)
  }
  proposed = proposed[proposed<upr & proposed>lwr]
  while(length(proposed)<samples){
    new_proposed=rnorm(samples-length(proposed))
    new_proposed = new_proposed[new_proposed<upr & new_proposed>lwr]
    proposed=c(proposed, new_proposed)
  }
  proposed[1:samples]
}, vectorize.args = c("upr", "lwr","loc","scale", "samples"))


generate_pd1_stacked = function(radiation_day, totaltime,pd1_times = c(-2,0,2,4)){
  num_pd1 = 4
  max_stack = 1.5
  
  pd1_mat = matrix(rep(0, (totaltime+7)*num_pd1), nrow = totaltime+7, ncol = num_pd1)
  for(pd1_idx in 1:num_pd1){
    start_time= pd1_times[pd1_idx]+15
    duration=7
    pdl1time = min(start_time+duration, totaltime)
    pd1_mat[(start_time):(pdl1time),pd1_idx] = pd1_mat[(start_time):(pdl1time),pd1_idx]+ seq(1,0, length.out = duration+1)[1:length(start_time:pdl1time)]
  }
  if(length(radiation_day)>0){
    for(pd1_idx in 1:num_pd1){
      start_times= pd1_times[pd1_idx]+radiation_day
      duration=7
      for(start_time in start_times){
        pdl1time = min(start_time+duration, totaltime)
        pd1_mat[(start_time):(pdl1time),pd1_idx] = pd1_mat[(start_time):(pdl1time),pd1_idx]+ seq(1,0, length.out = duration+1)[1:length(start_time:pdl1time)]
      }
    }
  }
  vec = as.numeric(rowSums(pd1_mat))
  vec[vec>max_stack]=max_stack
  return(vec)
}

SnT = function(d, alpha, beta, gam, Rn){
  logsnt = (1+gam*Rn)*(-alpha*d -beta*(d^2))
  return(exp(logsnt))
}

Rn_update = function(tau, SnT, Rn){
  vec = c(Rn-Rn/tau + (1-SnT),1)
  return(min(vec))
}

BTn_update= function(lambda, BTn, phi, d, rho){
  (1-lambda)*BTn*exp(-phi*d)+ rho*sign(d)
}

TTn_update = function(lambda, BTn, TTn, alpha1, d){
  lambda*BTn + TTn*exp(-alpha1*d)
}

Zn = function(omega1, omega2, TTn, BTn, p1){
  (p1)*omega1*TTn + (p1)*BTn*omega2
}

Tn_update = function(SnT, Tn, mu, Zn){
  SnT*Tn*exp(mu-Zn)
}
loglinear.growth=Vectorize(function(t, beta0, beta1){
  exp(beta0+beta1*t)
}, vectorize.args = "t")

tumor.decay = Vectorize(function(tt, dd, gg){
  exp(-gg*(tt-dd))
}, vectorize.args="tt")

library(boot)
get.mean.sumexp = Vectorize(function(rho1,theta1, beta0,beta1,gg, d1,current.time){
  scale.fact=1
  d0=sqrt(15)/scale.fact
  d1 = sqrt(d1)/scale.fact
  rhoo = rho1
  beta1 = beta1
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
}, vectorize.args = c("rho1","theta1","beta0","beta1","gg", "d1","current.time"))

generate_one_sumexp = function(radiation_days, parameter_vec, maxtime, add.measurement.noise=F){
  radiation_days=radiation_days[radiation_days>15]
  scale.fact=1
  d0=sqrt(15)/scale.fact
  if(length(radiation_days)>0){
    d1 = sqrt(radiation_days[1])/scale.fact
    }else{d1=sqrt(maxtime+1)/scale.fact}
  rhoo = parameter_vec['rho1']
  #*(1-min(c(exp(-parameter_vec['theta1']*(d0)-(0.3*parameter_vec['beta0']*d0)^2 ),1)))
  beta0 = parameter_vec['beta0']
  beta1 = parameter_vec['beta1'] 
  gg = abs(parameter_vec['gg'])
  rho2=parameter_vec['rho1']*(1-exp(-abs(parameter_vec['theta1'])*(d1-d0)))
  #print(rho2)
  time.vec=sqrt(1:maxtime)/scale.fact
  first.decay.term = rhoo*exp(beta0+beta1*d0)*exp(-gg*(time.vec-d0))*(time.vec>d0)
  fct.mean.exp = first.decay.term + (1-rhoo*(time.vec>d0))*exp(beta0+beta1*time.vec)
  if(maxtime>d1){
    fct.mean.exp.1=fct.mean.exp
    fct.mean.exp = rho2*(time.vec>d1)*fct.mean.exp.1*exp(-gg*(time.vec-d1))+(1-rho2*(time.vec>d1))*fct.mean.exp.1
    
  }else{
    
  }
  #log(fct.mean.exp)
  dose.vec = rep(0, length(fct.mean.exp))
  dose.vec[unique(c(15, radiation_days[radiation_days<maxtime]))] = 1
  out.array =array(data=NA, dim=c(1, length(fct.mean.exp),4))
  dimnames(out.array) = list(NULL, NULL, c("ltv","d","p",'day'))
  
  if(add.measurement.noise){
    noise.vec=truncnorm(1, lwr = -fct.mean.exp+.01, scale=00)
    fct.mean.exp =fct.mean.exp
  }
  out.array[1,,'ltv'] = log(fct.mean.exp)
  if(is.nan(sum(log(fct.mean.exp)))){
    print(parameter_vec)
    print(d1)
  }
  #print(dim(out.array))
  #print(length(cumsum(dose.vec)))
  out.array[1,,'d'] =cumsum(dose.vec)
  out.array[1,,'p'] = 0
  out.array[1,,'day'] = 1:maxtime
  out.array 
}


generate_one = function(radiation_days, parameter_vec, maxtime, gen.mod, add.measurement.noise=F){
  
  if(gen.mod =="recursive"){
  #print(rnorm(1, mean=0,sd=sqrt(init_variance)))
  para_set = make_default_para_set()
  para_set[names(parameter_vec),"best"] = as.numeric(parameter_vec)
  Tn_vec = rep(NA, maxtime)
  ###set vectors for model components
  SnT_vec<-rep(NA, maxtime)
  Zn_vec<-rep(NA, maxtime)
  TTn_vec=rep(NA, maxtime)
  BTn_vec=rep(NA, maxtime)
  Rn_vec<-rep(NA,maxtime)
  Rn_vec[1]<-0
  TTn_vec[1]<0
  TTn_vec[1]<-para_set['TT0','best']
  BTn_vec[1]<-para_set['BT0','best']
  dose_vec=rep(0, maxtime)
  dose_vec[15] = 10
  dose_vec[radiation_days[radiation_days>0]] = 10
  #print(radiation_days)
  pd1_vec = generate_pd1_stacked(radiation_days[radiation_days>0], totaltime=maxtime)
  
  
  Tn_vec[1]<- para_set["Tinit","best"]
  
  
  #print(Tn_vec[1])
  for (i in 1: (maxtime-1))
  { #print("looping...")
    SnT_vec[i]<-SnT(d=dose_vec[i], alpha=para_set['alpha','best'], beta=para_set['beta','best'],gam=para_set['gam','best'],Rn=Rn_vec[i])
    Rn_vec[i+1]<-Rn_update(tau=para_set['tau','best'],SnT=SnT_vec[i],Rn=Rn_vec[i])
    BTn_vec[i+1]<-BTn_update(lambda=para_set['lambda','best'], BTn=BTn_vec[i], phi=para_set['phi','best'], d=dose_vec[i], rho=para_set['rho','best'])
    TTn_vec[i+1]<-TTn_update(lambda=para_set['lambda','best'], BTn=BTn_vec[i], TTn=TTn_vec[i], alpha1=para_set['alpha1','best'],d=dose_vec[i])
    Zn_vec[i]<-Zn(omega1=para_set['omega1','best'], omega2=para_set['omega2','best'], TTn=TTn_vec[i], BTn=BTn_vec[i], p1=pd1_vec[i])
    Tn_vec[i+1]<-Tn_update(SnT=SnT_vec[i],Tn=Tn_vec[i],mu=para_set['mu','best'],Zn=Zn_vec[i])
  }
  #print(Tn_vec)
  #print(Tn_vec)
  #out.list = list(log(Tn_vec), dose_vec, pd1_vec)
  out.array =array(data=NA, dim=c(1, length(Tn_vec),4))
  dimnames(out.array) = list(NULL, NULL, c("ltv","d","p",'day'))
  out.array[1,,'ltv'] = log(Tn_vec)
  out.array[1,,'d'] =cumsum(dose_vec/10)
  out.array[1,,'p'] = pd1_vec[1:length(Tn_vec)]
  out.array[1,,'day'] = 1:length(Tn_vec)}else{
    
    out.array = generate_one_sumexp(radiation_days, parameter_vec, maxtime, add.measurement.noise=add.measurement.noise)
  }
  #names(out.list) = c("ltv","d","p")
  return(out.array)
}

make_parameter_mat = function(num_mice, gen.mod){
  if(gen.mod=="recursive"){
  para_set = make_default_para_set()
  mu=0.21601
  rho = 1.707
  lambda = 0.30441
  omega1 = 0.02349
  omega2 = 0.01404
  BT0 = 0.00001
  TT0 = 0
  Tinit = 1
  parameter_vec=c(mu, rho, lambda, omega1, omega2, BT0, TT0, Tinit)
  parameter_mat = matrix(rep(parameter_vec, num_mice), nrow=num_mice, byrow = T)
  colnames(parameter_mat) = c("mu","rho","lambda","omega1","omega2","BT0","TT0","Tinit")
  one.loc = lambda
  #parameter_mat[,'lambda'] = rbeta(num_mice,shape1 = one.loc*one.scale, shape2 = (1-one.loc)*one.scale)
  #parameter_mat[,'rho'] = sample(c(0,rho), num_mice, replace=T)
  parameter_mat[,"lambda"] = truncnorm(num_mice, loc=lambda,scale=.3, lwr=0, upr=1)
  #parameter_mat[,'BT0'] = truncnorm(num_mice, loc=BT0,scale=.001, lwr=0, upr=1)
  #parameter_mat[,"lambda"] = log(seq(1,exp(1), length.out=num_mice))
  parameter_mat[,"rho"] = truncnorm(num_mice, loc = rho, scale=.01,upr=3, lwr=0)
  parameter_mat[,"mu"] = truncnorm(num_mice, loc = mu, scale=.003,upr=1, lwr=0)
  #parameter_mat[,"omega1"] = truncnorm(num_mice, loc = omega1, scale=4e-6,upr=1, lwr=0.001)
  other_params = setdiff(rownames(para_set), colnames(parameter_mat))
  other_param_mat = matrix(rep(para_set[other_params,'best'], num_mice), ncol=length(other_params), byrow=T)
  colnames(other_param_mat)=other_params
  out.param.mat = cbind(parameter_mat, other_param_mat)}else{
    param_names = param_names=c("rho1", "theta1", "beta0","beta1","gg")
    minibatch_parameters = matrix(nrow=test_num, ncol=length(param_names))
    colnames(minibatch_parameters) = param_names
    true.param.gen = matrix(nrow=length(param_names), ncol=3)
    colnames(true.param.gen) = c("m","var",'upper')
    rownames(true.param.gen)=param_names
    true.param.gen[,'m'] = c(.9,.9,.1,.2,.2)
    true.param.gen[,'var'] =c(0.00,0,0,0.0000,1)
    true.param.gen[,'upper'] = c(1,1,1,1,1)
    for(pp in param_names){
      minibatch_parameters[,pp] = truncnorm(test_num, loc = true.param.gen[pp,'m'], scale=true.param.gen[pp,'var'], lwr=0, upr=true.param.gen[pp,'upper'])
    }
    out.param.mat = minibatch_parameters
    #out.param.mat[,'theta1'] = logit(out.param.mat[,'theta1'])
  }
  
  out.param.mat
}

generate_one_counterfactualset = function(parameter_vec, total_days, potential_actions, gen.mod){
  counterfactual_rtdays = 15+potential_actions
  counterfactualset.pairs = array(dim=c(length(potential_actions), total_days,4))
  
  dimnames(counterfactualset.pairs) = list(NULL, NULL, c("ltv","d","p",'day'))
  for(one.counterfactual.idx in 1:length(potential_actions)){
    counterfactualset.pairs[one.counterfactual.idx,,] = generate_one(radiation_days = c(counterfactual_rtdays[one.counterfactual.idx]),parameter_vec = parameter_vec, maxtime = total_days, gen.mod=gen.mod)
  
    }
  counterfactualset.pairs
}
