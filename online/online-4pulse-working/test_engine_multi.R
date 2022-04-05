setwd("~/RPAI/online/online-3pulse-working")
source("data_generation_fncs.R")
source("env_fncs.R")
setwd("~/RPAI/online/online-4pulse-working")

set.seed(1998)
test_num=50
maxtime = 120
num_free_pulses=3
wait_time=9
potential_actions = 1:13 + wait_time
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
inc_days = 15-2+wait_time
parameter_mat=make_parameter_mat(test_num)


test_num = 50

selected_actions=c()

agent.action.mat = matrix(NA, nrow=test_num, ncol=num_free_pulses)


agent.outcome.vec=rep(NA, test_num)
refdays = c(10,14,20)
ref.outcome.mat = matrix(NA, nrow=test_num, ncol=length(refdays)+1)
references=c(paste("day", refdays), "random")
colnames(ref.outcome.mat) = c(paste("day", refdays), "random")
est.par.mat = matrix(NA, nrow=test_num, ncol=3)
colnames(est.par.mat) = c("mu","lambda","rho")

for(mouse in 1:test_num){
  cat("mouse: ", mouse, "\n")
  parameter_vec=parameter_mat[mouse,]
  action.vec=c()
  radiation_days=c()
  one.sequence = generate_one(action.vec, parameter_vec, maxtime=20)$ltv
  while(length(action.vec)<num_free_pulses){
    isdone=(length(action.vec)==(num_free_pulses-1))
    estimated.par = optimize.model(one.sequence = one.sequence, radiation_days = radiation_days)
    names(estimated.par) = c("mu","lambda","rho")
    est.par.mat[mouse,] = estimated.par
    pulsno = length(action.vec)
    if(pulsno>0){
      
      selected.rows = 1:nrow(all.action.mat)
      for(jj in 1:length(action.vec)){
        selected.rows = intersect(selected.rows, which((all.action.mat[,jj]==action.vec[jj])))
      }
      rest.action.mat = all.action.mat[selected.rows,]
    }else{
      rest.action.mat=all.action.mat
    }
    one.plan = get.optim.plan(parameter_vec=estimated.par,maxtime=maxtime, all.action.mat = rest.action.mat)
    one.action=as.numeric(one.plan)[pulsno+1]
    action.vec=c(action.vec, one.action)
    radiation_days = cumsum(action.vec)+15
    seqlength = (inc_days+sum(action.vec))*(1-isdone)+maxtime*isdone
    one.sequence = generate_one(cumsum(action.vec)+15, parameter_vec, maxtime = seqlength)$ltv
  }
  for(ref.idx in 1:length(refdays)){
    ref = refdays[ref.idx]
    one.reference = generate_one(cumsum(rep(ref, num_free_pulses))+15, parameter_vec, maxtime = seqlength)$ltv
    ref.outcome.mat[mouse,ref.idx]=log(one.reference[seqlength])
  }
  one.random = sample(potential_actions,num_free_pulses, replace=T)
  ref.outcome.mat[mouse,length(refdays)+1]= log(generate_one(cumsum(one.random)+15, parameter_vec, maxtime = seqlength)$ltv)[seqlength]
  agent.action.mat[mouse,]=action.vec
  agent.outcome.vec[mouse] = log(one.sequence[seqlength])
  
}
effect.sizes = rep(NA, ncol(ref.outcome.mat))
names(effect.sizes) = references


x11()
par(mfrow=c(2,2))
adj.val=1
plot(density(-agent.outcome.vec[1:test_num]+ref.outcome.mat[1:test_num,1], adjust=adj.val), col="red", xlab="final ltv reference-final ltv agent",main="4 pulse performance")
colors= c("red","blue", "orange", "purple")
for(ref.idx in 1:ncol(ref.outcome.mat)){
  ref = refdays[ref.idx]
  lines(density(ref.outcome.mat[1:test_num,ref.idx] - agent.outcome.vec[1:test_num], adjust=adj.val), col=colors[ref.idx])
  print(t.test(ref.outcome.mat[,ref.idx] - agent.outcome.vec))
  delta = ref.outcome.mat[,ref.idx] - agent.outcome.vec
  delta2 = -ref.outcome.mat[,ref.idx] + ref.outcome.mat[,'random']
  effect.sizes[[references[ref.idx]]] = mean(delta)/sd(delta)
  print(summary(delta))
  print(mean(delta[delta<=0]))
  print(mean(delta[delta>0]))
  cat("proportion of times agent was better than", colnames(ref.outcome.mat)[ref.idx],"day spacing:" )
  print(1-length(delta[delta<=0])/100)
  cat("proportion of times ", colnames(ref.outcome.mat)[ref.idx],"was better than random" )
  print(1-length(delta[delta2<=0])/100)
}
legend("topright", legend=paste(c(paste("day", refdays),"random"), round(effect.sizes, digits=3), sep=": d="), pch=16, col=colors)
abline(v=0, lty=2)
mean(-agent.outcome.vec+ref.outcome.mat[,1])

plot(parameter_mat[,'lambda'], est.par.mat[,'lambda'],pch=16, main="lambda", ylab="estimated", xlab="actual")
plot(parameter_mat[,'rho'], est.par.mat[,'rho'], pch=16,main="rho", ylab="estimated", xlab="actual")
plot(parameter_mat[,'mu'], est.par.mat[,'mu'], pch=16,main="mu", ylab="estimated", xlab="actual")

