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
truncnorm = function(samples, loc=0, scale=1, upr=Inf, lwr=-Inf){
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
}


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


generate_one = function(radiation_days, parameter_vec, maxtime){
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
  dose_vec[radiation_days] = 10
  #print(radiation_days)
  pd1_vec = generate_pd1_stacked(radiation_days, totaltime=maxtime)
  
  
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
  out.list = list(log(Tn_vec), dose_vec, pd1_vec)
  out.array =array(data=NA, dim=c(1, length(Tn_vec),4))
  dimnames(out.array) = list(NULL, NULL, c("ltv","d","p",'day'))
  out.array[1,,'ltv'] = log(Tn_vec)
  out.array[1,,'d'] =cumsum(dose_vec/10)
  out.array[1,,'p'] = pd1_vec[1:length(Tn_vec)]
  out.array[1,,'day'] = 1:length(Tn_vec)
  names(out.list) = c("ltv","d","p")
  return(out.array)
}

make_parameter_mat = function(num_mice){
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

  parameter_mat[,"lambda"] = truncnorm(num_mice, loc=lambda,scale=.3, lwr=0, upr=1)
  #parameter_mat[,"lambda"] = rexp(num_mice, rate=1/lambda)

  #parameter_mat[,"rho"] = rexp(num_mice, rate = 1/rho)
  parameter_mat[,"rho"] = truncnorm(num_mice, loc = rho, scale=.01,upr=3, lwr=0)
  #parameter_mat[,"mu"] = rexp(num_mice, loc = mu, scale=.000025,upr=1, lwr=0)
  parameter_mat[,"omega1"] = truncnorm(num_mice, loc = omega1, scale=4e-6,upr=1, lwr=0.001)
  parameter_mat
}


