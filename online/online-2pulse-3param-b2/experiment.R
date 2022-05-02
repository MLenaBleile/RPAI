setwd("~/RPAI/online/online-2pulse-3param-b2")
source("data_generation_fncs.R")
source("env_fncs.R")
source("lm_fncs.R")


set.seed(1998)
num_free_pulses=1
total_time = 40
wait_time=9
potential_actions = 1:5+wait_time
train_idx=1:5
adapt_idx=1:85 + max(train_idx)
nonadapt_idx=0 + max(adapt_idx)
test_idx = (max(nonadapt_idx)-50):(max(nonadapt_idx))
total_mice=max(nonadapt_idx)
num_pc = 4
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)

##calibrate on initial data
parameter_mat = make_parameter_mat(total_mice)
action_mat = matrix(sample(potential_actions, total_mice, replace=T),byrow=T, nrow=total_mice, ncol=num_free_pulses)
action_mat[1:nrow(all.action.mat)] = all.action.mat[1:nrow(all.action.mat),]
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)


##do adaptive part
optim.acts= c()
rteffects = matrix(nrow=total_mice, ncol=7)
d2.vec = c()
for(one.mouse in 1:total_mice){
  op.at=get.optim.plan(parameter_mat[one.mouse,],40, all.action.mat)
  optim.acts = c(optim.acts, op.at)
  rteffects[one.mouse,] = get_rteffects(one_batch[one.mouse,1:22,])
  d2.vec = c(d2.vec, min(which(one_batch[one.mouse,,'d']==2)))
}
act.vec = d2.vec-15
for(nn in 1:7){
plot(rteffects[,nn]~optim.acts, pch=16, main=nn)
}

#fit = lm(optim.acts~(.)^2, data=as.data.frame(rteffects))
#table(round(predict(fit)),optim.acts)

pca.data = as.data.frame(rteffects)[train_idx,]
pca.obj = stats::prcomp(pca.data, scale=T)
in.data = as.data.frame(pca.obj$x[,1:num_pc])
#in.data$act = act.vec-9
in.data$act = as.character(act.vec)[train_idx]
actions = act.vec[train_idx]
fit2 = lm(one_batch[train_idx,40,'ltv']~(.)^2, data=in.data[train_idx,])

predicted.best = c()


for(one.mouse in adapt_idx){
  
  one.effects = get_rteffects(one_batch[one.mouse,1:22,])
  one.in = matrix(rep(as.numeric(one.effects), length(potential_actions)), byrow=T, ncol=7)
  one.in.df = as.data.frame(predict(pca.obj, as.data.frame(one.in))[,1:num_pc])
  one.in.df$act=as.character(potential_actions)

  one.best.pred = potential_actions[which.min(predict(fit2, one.in.df))]
  predicted.best = c(predicted.best, one.best.pred)
  pca.data = as.data.frame(rteffects)[1:one.mouse,]
  pca.obj = stats::prcomp(pca.data, scale=T)
  in.data = as.data.frame(pca.obj$x[,1:num_pc])
  actions = c(actions, one.best.pred)
  
  in.data$act = as.character(actions)
  one.paraset=matrix(parameter_mat[one.mouse,], nrow=1)
  colnames(one.paraset) = colnames(parameter_mat)
  one_batch[one.mouse,,] = generate_one_batch(parameter_mat = one.paraset, matrix(one.best.pred), current.time = total_time)[1,,]
  fit1 = glm(as.factor(act)~(.), family="binomial", data=in.data)
  data.weights = abs(predict(fit1) - actions)
  data.probs = data.weights/sum(data.weights)
  fit2 = lm(one_batch[1:one.mouse,40,'ltv']~(.)^2, data=in.data, weights=data.weights)

  
}



# for(one.mouse in nonadapt_idx){
#   #print("a")
#   one.effects = get_rteffects(one_batch[one.mouse,1:22,])
#   one.in = matrix(rep(as.numeric(one.effects), length(potential_actions)), byrow=T, ncol=7)
#   one.in.df = as.data.frame(one.in)
#   one.in.df$act=as.character(potential_actions)
#   
#   one.best.pred = potential_actions[which.min(predict(fit2, one.in.df[,1:8]))]
#   predicted.best = c(predicted.best, one.best.pred)
#   in.data$act[one.mouse] = one.best.pred
#   one.paraset=matrix(parameter_mat[one.mouse,], nrow=1)
#   colnames(one.paraset) = colnames(parameter_mat)
#   one_batch[one.mouse,,] = generate_one_batch(parameter_mat = one.paraset, matrix(one.best.pred), current.time = total_time)[1,,]
#   #fit2 = lm(one_batch[1:one.mouse,40,'ltv']~(.)^2, data=in.data[1:one.mouse,4:8])
# }


print(table(predicted.best[test_idx-max(train_idx)], optim.acts[test_idx]))
source("test_engine.R")
