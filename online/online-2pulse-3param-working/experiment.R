all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
burn_in=500
parameter_mat = make_parameter_mat(burn_in)
action_mat = matrix(sample(potential_actions, burn_in, replace=T),byrow=T, nrow=burn_in, ncol=num_free_pulses)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)

optim.acts= c()
rteffects = matrix(nrow=burn_in, ncol=7)
d2.vec = c()
for(jj in 1:500 ){
  op.at=get.optim.plan(parameter_mat[jj,],40, all.action.mat)
  optim.acts = c(optim.acts, op.at)
  rteffects[jj,] = get_rteffects(one_batch[jj,1:22,])
  d2.vec = c(d2.vec, min(which(one_batch[jj,,'d']==2)))
}
act.vec = d2.vec-15
for(nn in 1:7){
plot(rteffects[,nn]~optim.acts, pch=16, main=nn)
}

fit = lm(optim.acts~(.)^2, data=as.data.frame(rteffects))
table(round(predict(fit)),optim.acts)

in.data = as.data.frame(rteffects)
#in.data$act = act.vec-9
in.data$act = as.character(act.vec)
train_idx=1:40
fit2 = lm(one_batch[train_idx,40,'ltv']~(.)^2, data=in.data[train_idx,4:8])

predicted.best = c()
test_idx=41:120
for(jj in test_idx){
  
  one.effects = get_rteffects(one_batch[jj,1:22,])
  one.in = matrix(rep(as.numeric(one.effects), length(potential_actions)), byrow=T, ncol=7)
  one.in.df = as.data.frame(one.in)
  one.in.df$act=as.character(potential_actions)

  one.best.pred = potential_actions[which.min(predict(fit2, one.in.df[,4:8]))]
  predicted.best = c(predicted.best, one.best.pred)
  in.data$act[jj] = one.best.pred
  one.paraset=matrix(parameter_mat[jj,], nrow=1)
  colnames(one.paraset) = colnames(parameter_mat)
  one_batch[jj,,] = generate_one_batch(parameter_mat = one.paraset, matrix(one.best.pred), current.time = total_time)[1,,]
  fit2 = lm(one_batch[1:jj,40,'ltv']~(.)^2, data=in.data[1:jj,4:8])

  
  }
table(tail(predicted.best,80), tail(optim.acts[test_idx],80))
