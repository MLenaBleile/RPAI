setwd("~/RPAI/online/online-3pulse-3param-working")
source("data_generation_fncs.R")
source("env_fncs.R")
source("lm_fncs.R")
library(glmnet)


set.seed(1998)
num_free_pulses=2
total_time = 60
wait_time=9
potential_actions = 1:5+wait_time
train_num=800
adapt_num=0
test_num = 50
train_idx=1:(train_num*num_free_pulses)
adapt_idx=1:(adapt_num*num_free_pulses) + max(train_idx)
nonadapt_idx=1:(test_num*num_free_pulses) + max(adapt_idx)
test_idx = tail(c(adapt_idx,nonadapt_idx),100)
total_mice=train_num+adapt_num+test_num
num_pc = 7
inc.time=wait_time+15-2
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)

##calibrate on initial data
parameter_mat = make_parameter_mat(total_mice)
action_mat = all.action.mat[sample(1:nrow(all.action.mat), total_mice, replace=T),]
action_mat[1:nrow(all.action.mat),] = all.action.mat[1:nrow(all.action.mat),]
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)


##do adaptive part
optim.acts= matrix(nrow=total_mice, ncol=num_free_pulses)

one.pair= one_batch[1,1:total_time,]
test_effects = get_rteffects(one.pair)
rteffects = array(0, dim=c(total_mice, length(test_effects), 2))
dimnames(rteffects) = list(paste("mouse",1:total_mice, sep="_"),names(test_effects), c("d2","d3"))

dn.mat = matrix(NA, nrow=total_mice, ncol=num_free_pulses)
for(one.mouse in 1:total_mice){
  op.at=get.optim.plan(parameter_mat[one.mouse,],total_time, as.matrix(all.action.mat))
  optim.acts[one.mouse,] = op.at
  current.time=inc.time
  
  for(pp in 1:num_free_pulses){
  current.time = inc.time
  if(pp>1){
    current.time = min(which(one_batch[one.mouse,,'d']==(pp))) + wait_time
  }
  if(one.mouse<=train_num){
  one.pair= one_batch[one.mouse,1:current.time,]
  one.effectset = get_rteffects(one.pair)
  rteffects[one.mouse,,pp] = one.effectset
  dn.mat[one.mouse,pp] = min(which(one_batch[one.mouse,,'d']==(pp+1)))
    }
  }
}

act.mat = dn.mat -15
act.mat[,2] = act.mat[,2]-act.mat[,1]




library(dplyr)
rteffects.df = as.data.frame(rbind(rteffects[,,1], rteffects[,,2]))

pca.data = as.data.frame(rteffects.df)[train_idx,]
pca.obj = stats::prcomp(pca.data, scale=T)
in.data = as.data.frame(pca.obj$x[,1:num_pc])
#in.data$act = act.vec-9
actions = as.numeric(t(act.mat[1:train_num,]))
in.data$act = as.character(actions)[train_idx]
in.data$pulse = invlogit((rep(1:2, train_num)+.5)/3)
ltv.vec=rep(one_batch[1:train_num,total_time,'ltv'], num_free_pulses)
x <- model.matrix(ltv.vec ~ (.)^2, #-1 means no intercept
                  data = in.data[train_idx,])
fit2 = lm(ltv.vec~(.)^3, data=in.data[train_idx,])
#fit2 = lm(ltv.vec~(PC1+PC2+PC3+PC4+PC5+PC6+PC7)*act*pulse, data=in.data[train_idx,])

predicted.best = c()

obs.idx=length(actions)
#actions = as.matrix(actions, ncol=1)
names(actions) = paste(rep(1:train_num, each=num_free_pulses),rep(c(1:num_free_pulses), train_num), sep="_")

one.mouse=train_num+1
pp=1

# for(one.mouse in 1:adapt_num+train_num){
#   action.vec= c()
#   
#   for(pp in 1:num_free_pulses){
#     current.time = inc.time
#     if(pp>1){
#       current.time = min(which(one_batch[one.mouse,,'d']==(pp))) + wait_time
#     }
#   obs.idx= obs.idx+1
#   one.effects = get_rteffects(one_batch[one.mouse,1:current.time,])
#   rteffects[one.mouse,,pp] = one.effects
#   one.in = matrix(rep(as.numeric(one.effects), length(potential_actions)), byrow=T, ncol=length(test_effects))
#   colnames(one.in) = colnames(pca.data)
#   one.in.df = as.data.frame(predict(pca.obj, as.data.frame(one.in))[,1:num_pc])
#   one.in.df$act=as.character(potential_actions)
#   one.in.df$pulse=invlogit((pp+.5)/3)
# 
#   one.best.pred = potential_actions[which.min(predict(fit2, one.in.df))]
#   action.vec = c(action.vec, one.best.pred)
#   
#   predicted.best = c(predicted.best, one.best.pred)
#   actions = c(actions, one.best.pred)
#   names(actions)[length(actions)] = paste(one.mouse, pp, sep="_")
#   
#   pca.data.baseline = rbind(rteffects[1:(one.mouse-1),,1], rteffects[1:(one.mouse-1),,2])
#   pca.data.idx = paste(rep(1:(one.mouse-1), 2), c(rep(1, one.mouse-1), rep(2, one.mouse-1)), sep="_")
#   pca.data.new = matrix(rteffects[one.mouse,,1], nrow=1)
#   pca.data.idx = c(pca.data.idx, paste(one.mouse, 1, sep="_"))
#   colnames(pca.data.new) = dimnames(rteffects)[[2]]
#   if(pp==2){
#     pca.data.new= rbind(pca.data.new, rteffects[one.mouse,,2])
#     pca.data.idx = c(pca.data.idx, paste(one.mouse, pp, sep="_"))
#   }
#   pca.data= as.data.frame(rbind(pca.data.baseline, pca.data.new))
#   pca.obj = stats::prcomp(pca.data, scale=T)
#   in.data = as.data.frame(pca.obj$x[,1:num_pc])
#   
# 
#   in.data$act = actions[pca.data.idx]
#   one.paraset=matrix(parameter_mat[one.mouse,], nrow=1)
#   colnames(one.paraset) = colnames(parameter_mat)
#   one_batch[one.mouse,,] = generate_one_batch(parameter_mat = one.paraset, action_mat = matrix(action.vec, nrow=1), current.time = total_time)[1,,]
#   act.cts = invlogit((in.data$act+.5)/(max(in.data$act)+1))
#   fit1 = lm(act.cts~(.)^2, data=in.data)
#   data.weights = abs(predict(fit1) - invlogit((in.data$act+.5)/(max(in.data$act)+1)))
#   in.data$act = as.character(in.data$act)
#   in.data$pulse = c(rep(1, one.mouse-1), rep(2, one.mouse-1), 1:2)[1:length(actions)]
#   in.data$pulse =invlogit((in.data$pulse+.5)/3)
#   
#   data.probs = data.weights/sum(data.weights)
#   ltv.vec = c(one_batch[1:(one.mouse-1),total_time,'ltv'], one_batch[1:(one.mouse-1),total_time,'ltv'])
#   ltv.vec = c(ltv.vec,rep(one_batch[one.mouse,total_time,'ltv'],pp))
#   fit2 = lm(ltv.vec~(.)^3, data=in.data, weights=data.probs)
#   #fit2 = lm(ltv.vec~(PC1+PC2+PC3+PC4+PC5+PC6+PC7)*act*pulse, data=in.data, weights=data.probs)
#   }
#   act.mat[one.mouse,] = action.vec
# }



for(one.mouse in 1:test_num + train_num+adapt_num){
  action.vec= c()

  for(pp in 1:num_free_pulses){
    current.time = inc.time
    if(pp>1){
      current.time = min(which(one_batch[one.mouse,,'d']==(pp))) + wait_time
    }
    obs.idx= obs.idx+1
    one.effects = get_rteffects(one_batch[one.mouse,1:current.time,])
    rteffects[one.mouse,,pp] = one.effects
    one.in = matrix(rep(as.numeric(one.effects), length(potential_actions)), byrow=T, ncol=length(test_effects))
    colnames(one.in) = colnames(pca.data)
    one.in.df = as.data.frame(predict(pca.obj, as.data.frame(one.in))[,1:num_pc])
    one.in.df$act=as.character(potential_actions)
    one.in.df$pulse=invlogit((pp+.5)/3)

    one.best.pred = potential_actions[which.min(predict(fit2, one.in.df))]
    action.vec = c(action.vec, one.best.pred)

    predicted.best = c(predicted.best, one.best.pred)
    actions = c(actions, one.best.pred)
    names(actions)[length(actions)] = paste(one.mouse, pp, sep="_")
    #fit2 = lm(ltv.vec~(PC1+PC2+PC3+PC4+PC5+PC6+PC7)*act*pulse, data=in.data, weights=data.probs)
  }
  act.mat[one.mouse,] = action.vec
}

#test_idx = adapt_idx
optim.acts.long = as.numeric(t(optim.acts))[test_idx]
print(table(tail(predicted.best,length(optim.acts.long)) ,optim.acts.long))
test_mice= tail(train_num+1:adapt_num, 50)
print(table(act.mat[test_mice,1], optim.acts[test_mice,1]))
print(table(act.mat[test_mice,2], optim.acts[test_mice,2]))
print(head(act.mat))
source("test_engine.R")
