
generate_one_batch = function(minibatch_size,minibatch_parameters,wait_time, animal.idx=0){
  maxtime=15-2+wait_time
  data_mat = matrix(NA,nrow=0, ncol=4)
  colnames(data_mat) = c("ltv", "day","animal", "action")
  for(mouse.idx  in 1:minibatch_size){
    one.sequence = generate_one(c(), parameter_vec = minibatch_parameters[mouse.idx,], maxtime=maxtime)$ltv
    one.datum = as.matrix(data.frame(ltv= log(one.sequence), day=1:maxtime, animal = animal.idx+mouse.idx, action=rep(0, maxtime)))
    data_mat = rbind(data_mat, one.datum)
  }
  data_mat
}

update_one_batch = function(data_mat, minibatch_parameters, action_mat, maxtime=NA){
  animal.indices = unique(data_mat[,'animal'])
  new_data_mat = matrix(NA,nrow=0, ncol=4)
  colnames(new_data_mat) = c("ltv", "day","animal", "action")
  for(mouse.idx in 1:minibatch_size){
    one.action.vec = action_mat[mouse.idx,]
    if(is.na(maxtime)){
      maxtime = 15-2+sum(one.action.vec)
    }
    radiation_days = cumsum(one.action.vec)+15 
    one.pair = generate_one(radiation_days, parameter_vec = minibatch_parameters[mouse.idx,], maxtime=maxtime)
    one.sequence=one.pair$ltv
    one.dose=one.pair$d
    action.factor = rep(0, maxtime)
    action.factor[cumsum(one.dose)==2] = rep(one.action.vec[1], sum(cumsum(one.dose)==2))
    one.datum = as.matrix(data.frame(ltv=log(one.sequence), 1:maxtime, animal= rep(animal.indices[mouse.idx], maxtime), action= action.factor))
    new_data_mat = rbind(new_data_mat, one.datum)
  }
  new_data_mat
}

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
obj.fnc = function(input_vec, pair.set,param_names, bounds, fixed_param_vec){
  num_sequences=dim(pair.set)[1]
  loss.contribs=numeric(num_sequences)
  parameter_vec = r.to.unit(x=input_vec)
  
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
    pd1_vec=one.pair[,'p']
    sim.mod = lm(one.sim~day.vec*dose.vec*pd1_vec)
    real.mod = lm(one.sequence~day.vec*dose.vec*pd1_vec)
    ss.sim = summary(sim.mod)
    ss.real = summary(real.mod)
    sim.info = c(ss.sim$r.squared, 2*as.numeric(ss.sim$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
    real.info = c(ss.real$r.squared, 2*as.numeric(ss.real$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
    #sqrt(mean((log(one.sim[mt])-log(one.sequence[mt]))^2))
    loss.contribs[sequence.idx] = sqrt(mean((sim.info-real.info)^2)) + (0)*sqrt(mean((one.sim - one.sequence)^2))
  }
  mean(loss.contribs)
}

optimize.model = function(pair.set, param_names, fixed_param_vec){
  input_vec = unit.to.r(runif(length(param_names)))
  names(input_vec) = param_names
  input_vec['rho'] = unit.to.r(runif(1, 0,3),3)
  bounds = matrix(nrow=7, ncol=2)
  colnames(bounds) = c('upper','lower')
  rownames(bounds) = c('omega1','omega2','alpha','beta','alpha1','phi','tau')
  bounds[,'lower'] = c(0.02,0.01,0,0,0,.5,5)
  bounds[,'upper'] = c(.03,.02,.1,.1,.1,1,10)
  for(rparam in intersect(rownames(bounds), param_names)){
    input_vec[rparam] = unit.to.r(runif(1, bounds[rparam,'lower'],bounds[rparam,'upper']),s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
  }
  optim.obj = optim(input_vec, obj.fnc, pair.set=pair.set, bounds=bounds, param_names=param_names,fixed_param_vec=fixed_param_vec)
  estimated.input = optim.obj$par
  estimated.par = r.to.unit(estimated.input)
  names(estimated.par) = names(estimated.input)
  estimated.par['rho'] = estimated.par['rho']*3
  for(rparam in intersect(rownames(bounds), param_names)){
    estimated.par[rparam] = r.to.unit(estimated.input[rparam],s=bounds[rparam,'upper'], shift=bounds[rparam, 'lower'])
  }
  estimated.par
}


