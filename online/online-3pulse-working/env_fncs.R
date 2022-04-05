get.optim.plan = function(parameter_vec,maxtime, all.action.mat){
  ftv.vec= rep(NA, nrow(all.action.mat))
  for(j in 1:nrow(all.action.mat)){
    one.action.vec = as.numeric(all.action.mat[j,])
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