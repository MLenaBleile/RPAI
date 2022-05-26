get.optim.plan = function(parameter_vec,maxtime, all.action.mat){
  ftv.vec= rep(NA, nrow(all.action.mat))
  for(j in 1:nrow(all.action.mat)){
    one.action.vec = all.action.mat[j,]
    one.sequence = generate_one(radiation_days=c(cumsum(one.action.vec)+15), parameter_vec=parameter_vec, maxtime=maxtime+10)[1,,'ltv']
    ftv.vec[j] = one.sequence[maxtime]
  }
  best.action.vec=all.action.mat[which.min(ftv.vec),]
  best.action.vec
}

r.to.unit = Vectorize(function(x, s=1, shift=0){(s-shift)/(1+exp(x))+shift}, vectorize.args = c("x","s"))
unit.to.r = Vectorize(function(y,s=1, shift=0){log((s-shift)/(y-shift) - 1)}, vectorize.args = c("y","s"))
obj.fnc = function(input_vec, pair.set,param_names, bounds, fixed_param_vec, liklihood.based){
  num_sequences=dim(pair.set)[1]
  loss.contribs=numeric(num_sequences)
  parameter_vec = r.to.unit(x=input_vec[param_names])
  if(num_sequences>1){
    sig.vals = r.to.unit(input_vec['sigma'])*sqrt(1:length(pair.set[1,,'ltv']))
    sig.vals[is.na(sig.vals)]=1
    #print(sig.vals)
  }else{sig.vals=rep(1, length(pair.set[1,,'ltv']))}
  
  names(parameter_vec) = param_names
  parameter_vec['rho'] = r.to.unit(input_vec['rho'])*3
  for(rparam in intersect(rownames(bounds), param_names)){
    parameter_vec[rparam] = r.to.unit(input_vec[rparam],s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
  }
  
  for(fixed_param_name in names(fixed_param_vec)){
    parameter_vec[fixed_param_name] = fixed_param_vec[fixed_param_name]
  }
  for(sequence.idx in 1:num_sequences){
    one.sequence = pair.set[sequence.idx,,'ltv']
    one.pair=pair.set[sequence.idx,,]
    mt=length(one.sequence)
    radiation_days = which(one.pair[,'d']==2)[1]
    one.sim.pair = generate_one(radiation_days = c(),parameter_vec=parameter_vec, maxtime = mt)
    one.sim = one.sim.pair[1,,'ltv']
    day.vec=one.pair[,'day']
    dose.vec=one.pair[,'d']
    #pd1_vec=one.pair[,'p']
    sim.mod = lm(one.sim~day.vec*dose.vec)
    real.mod = lm(one.sequence~day.vec*dose.vec)
    ss.sim = summary(sim.mod)
    ss.real = summary(real.mod)
    sim.info = c(ss.sim$r.squared, 2*as.numeric(ss.sim$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
    real.info = c(ss.real$r.squared, 2*as.numeric(ss.real$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
    #sqrt(mean((log(one.sim[mt])-log(one.sequence[mt]))^2))
    one.standardized = one.sequence - one.sim
    lik.contrib = log(dnorm(one.standardized, sd=sig.vals)+.001)
    if(is.infinite(sum(lik.contrib))){cat(dnorm(one.standardized, sd=sig.vals))}
    loss.contribs[sequence.idx] = (1-liklihood.based)*sqrt(mean((sim.info-real.info)^2)) - sum(liklihood.based*lik.contrib)
  }
  mean(loss.contribs)
}

optimize.model = function(pair.set, param_names, fixed_param_vec, starting.params=NA,liklihood.based=F){
  input_vec = unit.to.r(runif(length(param_names)))
  names(input_vec) = param_names
  input_vec['rho'] = unit.to.r(runif(1, 0,3),3)
  bounds = matrix(nrow=7, ncol=2)
  colnames(bounds) = c('upper','lower')
  rownames(bounds) = c('omega1','omega2','alpha','beta','alpha1','phi','tau')
  bounds[,'lower'] = c(0.02,0.01,0,0,0,.5,5)
  bounds[,'upper'] = c(.03,.02,.1,.1,.1,1,10)
  if(liklihood.based){
    #print("doing liklihood based calculation")
    #input_vec = c(input_vec, unit.to.r(runif(1)))
    input_vec = starting.params
    if(sum(is.na(input_vec))>0){
      print("randomly initializing sigma")
      input_vec[is.na(input_vec)] = runif(1)
    }
    names(input_vec) = c(param_names, "sigma")
  }
  for(rparam in intersect(rownames(bounds), param_names)){
    input_vec[rparam] = unit.to.r(runif(1, bounds[rparam,'lower'],bounds[rparam,'upper']),s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
  }
  #print(liklihood.based)
  optim.obj = optim(input_vec, obj.fnc, pair.set=pair.set, bounds=bounds, param_names=param_names,fixed_param_vec=fixed_param_vec, liklihood.based=liklihood.based)
  estimated.input = optim.obj$par
  estimated.par = r.to.unit(estimated.input)
  names(estimated.par) = names(estimated.input)
  estimated.par['rho'] = estimated.par['rho']*3
  for(rparam in intersect(rownames(bounds), param_names)){
    estimated.par[rparam] = r.to.unit(estimated.input[rparam],s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
  }
  estimated.par
}


