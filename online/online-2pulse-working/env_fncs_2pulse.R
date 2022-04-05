
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
    one.sequence = generate_one(radiation_days=c(cumsum(one.action.vec)+15), parameter_vec=parameter_vec, maxtime=maxtime+10)$ltv
    ftv.vec[j] = one.sequence[maxtime]
  }
  best.action.vec=all.action.mat[which.min(ftv.vec),]
  best.action.vec
}

r.to.unit = function(x, s=1){s/(1+exp(x))}
unit.to.r = function(y,s=1){log(s/y - 1)}
obj.fnc = function(input_vec, one.sequence,radiation_days=c()){
  mt=length(one.sequence)
 
  parameter_vec = c(r.to.unit(input_vec[1]), r.to.unit(input_vec[2]), (r.to.unit(input_vec[3])*.5+1.5))
  names(parameter_vec) = c("mu","lambda","rho")
  one.sim.pair = generate_one(radiation_days = radiation_days,parameter_vec=parameter_vec, maxtime = mt)
  one.sim = one.sim.pair$ltv
  day.vec=1:mt
  dose.vec=cumsum(one.sim.pair$d/10)
  pd1_vec=one.sim.pair$p[1:length(day.vec)]
  sim.mod = lm(log(one.sim)~day.vec*dose.vec*pd1_vec)
  real.mod = lm(log(one.sequence)~day.vec*dose.vec*pd1_vec)
  ss.sim = summary(sim.mod)
  ss.real = summary(real.mod)
  sim.info = c(ss.sim$r.squared, 2*as.numeric(ss.sim$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
  real.info = c(ss.real$r.squared, 2*as.numeric(ss.real$coefficients[,1]), as.numeric(ss.real$coefficients[,2]))
  #sqrt(mean((log(one.sim[mt])-log(one.sequence[mt]))^2))
  sqrt(mean((sim.info-real.info)^2)) + (0)*sqrt(mean((log(one.sim) - log(one.sequence))^2))
}

optimize.model = function(one.sequence, radiation_days=c()){
  input_vec = c(unit.to.r(runif(1)),unit.to.r(runif(1)), unit.to.r(runif(1,0,3),3))
  optim.obj = optim(input_vec, obj.fnc, one.sequence=one.sequence, radiation_days=radiation_days)
  estimated.input = optim.obj$par
  estimated.par = c(r.to.unit(estimated.input[1]), r.to.unit(estimated.input[2]), r.to.unit((estimated.input[3]))*.5+1.5)
  estimated.par
}


