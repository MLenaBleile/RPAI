setwd("~/RPAI/online/online-general-lm - recursive")
source("data_generation_fncs.R")
source("env_fncs.R")
source("lm_fncs.R")
library(glmnet)


set.seed(1998)
num_free_pulses=1
total_time = 40 + 20*(num_free_pulses-1)
wait_time=9
potential_actions = 1:5+wait_time
train_num=10/num_free_pulses
adapt_num=100/num_free_pulses
test_num = 50
train_idx=1:(train_num*num_free_pulses)
adapt_idx=1:(adapt_num*num_free_pulses) + max(train_idx)
nonadapt_idx=1:(test_num*num_free_pulses) + max(adapt_idx)
test_idx = tail(c(adapt_idx),test_num*num_free_pulses)
total_mice=train_num+adapt_num
num_pc = 7
inc.time=wait_time+15-2
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)

##make initial data
parameter_mat = make_parameter_mat(total_mice)
action_mat = all.action.mat[sample(1:nrow(all.action.mat), total_mice, replace=T),]
if(num_free_pulses==1){
  action_mat = matrix(action_mat, ncol=1)
}
action_mat[1:nrow(all.action.mat),] = as.matrix(all.action.mat)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)


##get optimal actions and features for training
optim.acts= matrix(nrow=total_mice, ncol=num_free_pulses)

one.pair= one_batch[1,1:total_time,]
test_effects = get_rteffects(one.pair, num_free_pulses, total_time)
rteffects = array(0, dim=c(total_mice, length(test_effects)+1, num_free_pulses))
dimnames(rteffects) = list(1:total_mice,c(names(test_effects), 'action'), as.character(1:num_free_pulses+1))
names(dimnames(rteffects)) = c("animal", "effect", "pulse")

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
  one.effectset = get_rteffects(one.pair, num_free_pulses = num_free_pulses, total_time = current.time)
  rteffects[one.mouse,1:length(test_effects),pp] = one.effectset
  dn.mat[one.mouse,pp] = min(which(one_batch[one.mouse,,'d']==(pp+1)))
  rteffects[one.mouse,'action',pp] = min(which(one_batch[one.mouse,,'d']==(pp+1)))- min(which(one_batch[one.mouse,,'d']==(pp)))
    }
  }
}


##calibrate on initial data
library(dplyr)
library(reshape2)
rteffects.molten = melt(rteffects, names(dimnames(rteffects)))
all.pca.data = dcast(rteffects.molten, animal+pulse~effect, value.var = 'value')[train_idx,]
pca.data=all.pca.data[,-c(1,ncol(all.pca.data), which(apply(all.pca.data, 2, sd)==0))]
pca.obj = stats::prcomp(pca.data, scale=T)
in.data = as.data.frame(pca.obj$x[,1:num_pc])
actions = as.numeric(t(action_mat[1:train_num,]))
in.data$act = as.character(all.pca.data[,'action'])[train_idx]
ltv.vec=rep(one_batch[1:train_num,total_time,'ltv'], num_free_pulses)
in.data$ltv=ltv.vec
fit2 = lm(ltv~(.)^3, data=in.data[train_idx,])
#fit2 = nnet::nnet(ltv~(.)^3, data=in.data[train_idx,], size=5, linout=T)

predicted.best = c()

obs.idx=sum(rteffects[,'action',]!=0)

one.mouse=train_num+1
pp=1
##initialize df for counterfactuals
all.cf.preds = data.frame()

for(one.mouse in 1:adapt_num+train_num){
  cat("optimizing mouse,",one.mouse," \n")
  one.paraset=matrix(parameter_mat[one.mouse,], nrow=1)
  colnames(one.paraset) = colnames(parameter_mat)
  ###prediction phase
  for(pp in 1:num_free_pulses){
    current.time = inc.time
    if(pp>1){
      current.time = min(which(one_batch[one.mouse,,'d']==(pp))) + wait_time
    }else{
      one_batch[one.mouse,1:current.time,] = generate_one(c(),parameter_vec=one.paraset, current.time)
    }
  obs.idx= obs.idx+1
  one.effects = get_rteffects(one_batch[one.mouse,1:current.time,],num_free_pulses = num_free_pulses, total_time=current.time)
  rteffects[one.mouse,,pp] = c(one.effects,NA)
  if(num_free_pulses>1){
  one.in.vec = c(pp+1, one.effects)
  names(one.in.vec)[1]="pulse"}else{one.in.vec=c(one.effects)}
  one.in = matrix(rep(one.in.vec, length(potential_actions)), byrow=T, nrow=length(potential_actions))
  colnames(one.in) =names(one.in.vec)
  one.in.df = as.data.frame(predict(pca.obj, as.data.frame(one.in))[,1:num_pc])
  one.in.df$act=as.character(potential_actions)
  predictions=predict(fit2, one.in.df)
  pred.df = one.in.df
  pred.df$ltv = predictions
  selected.idx=which.min(predictions)
  one.best.pred = potential_actions[selected.idx]
  
  rteffects[one.mouse,'action',as.character(pp+1)] = one.best.pred
  

  ###Model refitting phase
  names.formelt=names(dimnames(rteffects))
  if(num_free_pulses==1){
    #names.formelt=names(dimnames(rteffects))
    names.formelt = names.formelt[1:2]
    rteffects.molten = melt(rteffects[1:(one.mouse),,], names.formelt)
    
    all.pca.data = dcast(rteffects.molten, animal~effect, value.var = 'value')
    mouse.idx = all.pca.data$animal
    pca.data= all.pca.data[,-c(1, ncol(all.pca.data), which(apply(all.pca.data, 2, sd)==0))]
  }else{
    #names.formelt=names(dimnames(rteffects))
    rteffects.molten.wextra = melt(rteffects[1:(one.mouse),,], names.formelt)
    #get rid of obs for current mouse that haven't happened yet
    if(pp<num_free_pulses){
      rteffects.molten = rteffects.molten.wextra[-c(which(rteffects.molten.wextra$pulse>(pp+1)&rteffects.molten.wextra$animal==one.mouse)),]}else{
        rteffects.molten=rteffects.molten.wextra
      }
    all.pca.data = dcast(rteffects.molten, animal+pulse~effect, value.var = 'value')
    mouse.idx = all.pca.data$animal
    pca.data= all.pca.data[,-c(1, ncol(all.pca.data), which(apply(all.pca.data, 2, sd)==0))]
  }
  
  
  pca.obj = stats::prcomp(pca.data, scale=T)
  in.data = as.data.frame(pca.obj$x[,1:num_pc])
  

  in.data$act = as.character(all.pca.data[,'action'])
  
  colnames(one.paraset) = colnames(parameter_mat)
  one_batch[one.mouse,,] = generate_one(parameter_vec = one.paraset, radiation_days = cumsum(rteffects[one.mouse,'action',])+15, maxtime=total_time)

  ltv.vec = rep(one_batch[1:(one.mouse-1),total_time,'ltv'], each=num_free_pulses)
  ltv.vec = c(ltv.vec,rep(one_batch[one.mouse,total_time,'ltv'],pp))
  in.data$ltv = ltv.vec
  
  all.cf.preds = rbind(all.cf.preds, pred.df[-c(selected.idx),])
  
  fit2 = lm(ltv~(.)^3, data=rbind(in.data, all.cf.preds))
  #fit2= nnet::nnet(ltv~(.)^3, data=rbind(in.data, all.cf.preds), size=5, linout=T)
  }

}



# for(one.mouse in 1:test_num + train_num+adapt_num){
#   action.vec= c()
#   
#   for(pp in 1:num_free_pulses){
#     current.time = inc.time
#     if(pp>1){
#       current.time = min(which(one_batch[one.mouse,,'d']==(pp))) + wait_time
#     }
#     obs.idx= obs.idx+1
#     one.effects = get_rteffects(one_batch[one.mouse,1:current.time,])
#     rteffects[one.mouse,,pp] = one.effects
#     one.in = matrix(rep(as.numeric(one.effects), length(potential_actions)), byrow=T, ncol=length(test_effects))
#     colnames(one.in) = colnames(pca.data)
#     one.in.df = as.data.frame(predict(pca.obj, as.data.frame(one.in))[,1:num_pc])
#     one.in.df$act=as.character(potential_actions)
#     one.in.df$pulse=invlogit((pp+.5)/3)
#     
#     one.best.pred = potential_actions[which.min(predict(fit2, one.in.df))]
#     action.vec = c(action.vec, one.best.pred)
#     
#     predicted.best = c(predicted.best, one.best.pred)
#     actions = c(actions, one.best.pred)
#     names(actions)[length(actions)] = paste(one.mouse, pp, sep="_")
#     #fit2 = lm(ltv.vec~(PC1+PC2+PC3+PC4+PC5+PC6+PC7)*act*pulse, data=in.data, weights=data.probs)
#   }
#   act.mat[one.mouse,] = action.vec
# }

source("test_engine.R")
