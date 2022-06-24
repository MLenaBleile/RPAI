true.param.gen = matrix(NA, nrow=length(param_names),ncol=3)
colnames(true.param.gen) = c("m","var",'upper')
rownames(true.param.gen) = param_names
true.param.gen[,'m'] = c(.9,.9,.1,.216,.5)
true.param.gen[,'var'] =c(.1,0,0.0,0.0,1)
true.param.gen[,'upper'] = c(1,1,1,1,1)
for(pp in param_names){
  minibatch_parameters[,pp] = truncnorm(test_num, loc = true.param.gen[pp,'m'], scale=true.param.gen[pp,'var'], lwr=0, upr=true.param.gen[pp,'upper'])
}

for(mouse in 1:test_num){
  cat("optimizing mouse ",mouse,"\n")
  true.param.vec= minibatch_parameters[mouse,]
  one.reference.pairs = generate_one_counterfactualset(minibatch_parameters[mouse,],total_days = maxtime,potential_actions = 1:20, gen.mod="sumexp")
  one.optima= potential_actions[which.min(one.reference.pairs[potential_actions,maxtime,'ltv'])]
  optimal_actions[mouse] = one.optima
}
print(table(optimal_actions))