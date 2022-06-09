
get.optim.plan = function(parameter_vec,maxtime, all.action.mat, mod.to.use){
  ftv.vec= rep(NA, nrow(all.action.mat))
  if(mod.to.use=="recursive"){
      for(j in 1:nrow(all.action.mat)){
        one.action.vec = all.action.mat[j,]
        one.sequence = generate_one(radiation_days=c(cumsum(one.action.vec)+15), parameter_vec=parameter_vec, maxtime=maxtime+10)[1,,'ltv']
        ftv.vec[j] = one.sequence[maxtime]
      }
  }else if(mod.to.use=="sumexp"){
      for(j in 1:nrow(all.action.mat)){
        one.action.vec = all.action.mat[j,]
        ftv.vec[j]=generate_one_sumexp(cumsum(one.action.vec)+15, parameter_vec=parameter_vec,maxtime=maxtime)[maxtime]
      }
    }
  best.action.vec=all.action.mat[which.min(ftv.vec),]
  best.action.vec
}

r.to.unit = Vectorize(function(x, s=1, shift=0){(s-shift)/(1+exp(x))+shift}, vectorize.args = c("x","s"))
unit.to.r = Vectorize(function(y,s=1, shift=0){log((s-shift)/(y-shift) - 1)}, vectorize.args = c("y","s"))
obj.fnc = function(input_vec, pair.set,param_names, bounds, fixed_param_vec, loss.weights, mod.to.fit, convert){
  
  num_sequences=dim(pair.set)[1]
  loss.contribs=numeric(num_sequences)
  #parameter_vec = r.to.unit(x=input_vec)
  random_names = setdiff(param_names, names(fixed_param_vec))
  if(length(random_names)==1){
    names(input_vec) = random_names}else{}
  parameter_vec = c(input_vec[random_names], fixed_param_vec)
  names(parameter_vec)[1:length(random_names)] = random_names
  #print(parameter_vec)
  if(loss.weights['liklihood']>0){
    sig.vals = r.to.unit(input_vec['sigma'],s=bounds["sigma","upper"]-bounds["sigma","lower"], shift=bounds["sigma","lower"])
  }else{sig.vals=rep(1, length(pair.set[1,,'ltv']))}
  
  #names(parameter_vec) = c(random_names, names(fixed_param_vec))
  #parameter_vec['rho'] = r.to.unit(input_vec['rho'])*3
  if(convert){
    for(rparam in intersect(rownames(bounds), param_names)){
      parameter_vec[rparam] = r.to.unit(input_vec[rparam],s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
    }
    for(fixed_param_name in names(fixed_param_vec)){
      parameter_vec[fixed_param_name] = fixed_param_vec[fixed_param_name]
    }
    
  }
  
  
  for(sequence.idx in 1:num_sequences){
    one.sequence = pair.set[sequence.idx,,'ltv']
    
    one.pair=pair.set[sequence.idx,!is.na(one.sequence),]
    one.sequence = one.sequence[!is.na(one.sequence)]
    mt=length(one.sequence)
    radiation_days = which(one.pair[,'d']==2)[1]
    if(is.na(radiation_days)){
      radiation_days=c(15)
    }

    if(mod.to.fit=="recursive"){one.sim.pair = generate_one(radiation_days = radiation_days,parameter_vec=parameter_vec, maxtime = mt)
    
    }else if(mod.to.fit=="sumexp"){
     
      one.sim.pair = generate_one_sumexp(radiation_days = radiation_days,parameter_vec=parameter_vec, maxtime = mt)
    }
    one.sim = one.sim.pair[1,,'ltv']
    #print(one.sim)
    day.vec=one.pair[1:length(one.sequence),'day']
    #print(day.vec)
    dose.vec=one.pair[,'d']
    pd1_vec=one.pair[,'p']
    #print(parameter_vec)

    sim.mod = lm(one.sim~day.vec*dose.vec)
    real.mod = lm(one.sequence~day.vec*dose.vec)
    ss.sim = summary(sim.mod)
    ss.real = summary(real.mod)
    sim.info = c(ss.sim$r.squared, 2*as.numeric(ss.sim$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
    real.info = c(ss.real$r.squared, 2*as.numeric(ss.real$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
    
    one.standardized = one.sequence - one.sim
    lik.contrib = -log(dnorm(one.standardized, sd=sig.vals))*loss.weights['liklihood']
    effect.contrib = loss.weights['effects']*sqrt(mean((sim.info-real.info)^2))
    ssq.contrib = loss.weights['raw']*sqrt(mean((one.sim - one.sequence)^2))
    #print(sum(lik.contrib))
    loss.contribs[sequence.idx] =  effect.contrib/num_sequences + ssq.contrib/num_sequences +sum(lik.contrib, na.rm=T)
  }
  all.loss=sum(loss.contribs, na.rm=T)
  all.loss
}

optimize.model = function(pair.set, param_names, fixed_param_vec, loss.weights, bounds, mod.to.fit){
  input_vec = unit.to.r(runif(length(param_names)))
  names(input_vec) = param_names
  convert=length(setdiff(param_names, names(fixed_param_vec)))>1
  if(loss.weights['liklihood']>0){
    input_vec = c(input_vec, unit.to.r(runif(1)))
    names(input_vec) = c(param_names, "sigma")
  }
  
  
  if(convert){
    for(rparam in intersect(rownames(bounds), param_names)){
      input_vec[rparam] = unit.to.r(runif(1, bounds[rparam,'lower'],bounds[rparam,'upper']),s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
    }
   
  optim.obj = optim(input_vec, obj.fnc, pair.set=pair.set, bounds=bounds,method="Nelder-Mead", mod.to.fit = mod.to.fit,loss.weights=loss.weights,param_names=param_names,fixed_param_vec=fixed_param_vec, convert=convert)
  estimated.input = optim.obj$par
  estimated.par = estimated.input
  names(estimated.par) = names(estimated.input)
  #estimated.par['rho'] = estimated.par['rho']*3
  for(rparam in intersect(rownames(bounds), param_names)){
    estimated.par[rparam] = r.to.unit(estimated.input[rparam],s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
  }
  }else{
    #print(input_vec)
    optim.obj = optimize(obj.fnc,interval=c(0,0.5), pair.set=pair.set, bounds=bounds, mod.to.fit = mod.to.fit,loss.weights=loss.weights,param_names=param_names,fixed_param_vec=fixed_param_vec, convert=convert)
    estimated.par = optim.obj$minimum
    names(estimated.par) = setdiff(param_names,names(fixed_param_vec))
    #names(estimated.par) = names(estimated.input)
    }
    
  
  
  estimated.par
}


