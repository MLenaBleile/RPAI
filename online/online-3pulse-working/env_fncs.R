make_potential_action_mat = function(potential_actions, num_free_pulses){
all.action.mat = expand.grid(rep(list(potential_actions), num_free_pulses))
rownames(all.action.mat) = get_labels(all.action.mat = all.action.mat)
all.action.mat
}


generate_one_initial_batch = function(parameter_mat,inc.time, total_time){
  minibatch_size = nrow(parameter_mat)
  one_batch = matrix(NA, nrow=minibatch_size, ncol=total_time)
  for(idx in 1:minibatch_size){
    one.pair = generate_one(c(), parameter_mat[idx,], inc.time)
    one_batch[idx,1:inc.time] = one.pair$ltv
  }
  one_batch
}

generate_one_batch = function(parameter_mat, action_mat, current.times, total_time){
  minibatch_size = nrow(parameter_mat)
  
  one_batch = matrix(NA, nrow=minibatch_size, ncol=total_time)
  for(idx in 1:minibatch_size){
    active.time=min(current.times[idx], total_time)
    one.pair = generate_one(cumsum(action_mat[idx,])+15, parameter_mat[idx,], total_time)
    one_batch[idx,1:active.time] = one.pair$ltv[1:active.time]
  }
  one_batch
}

get_labels = function(all.action.mat){
  char.mat = apply(all.action.mat,c(2),as.character)
  char.labs = char.mat[,1]
  if(ncol(char.mat)>1){
  for(jj in 2:ncol(char.mat)){
    char.labs = paste(char.labs, char.mat[,jj], sep=".")
  }}
  char.labs
}

update_model = function(q.fit,all_data, one_batch, action_mat){
  character.actions = get_labels(action_mat)
  unique.trt = unique(character.actions)
  for(one.group.idx in 1:length(unique.trt)){
  one.group.mice = which(character.actions==unique.trt[one.group.idx])
  one_ds = as.data.frame(rbind(all_data[one.group.mice,],one_batch))
  long_days = 1:ncol(one_ds)
  vnames = colnames(one_ds)
  one_ds$animal = 1:nrow(one_ds)
  one_ds.long = reshape2::melt(one_ds, id.vars=c("animal"))
  one_ds.long$day = rep(long_days, length(one.group.mice)+nrow(one_batch))
  one.model = glm(value~as.factor(day), data=one_ds.long)
  q.fit[[unique.trt[one.group.idx]]] = one.model
  }
  q.fit
}

get_predictions = function(q.fit, one_batch, all.action.mat,action_mat,pulse, add_variability=T){
  prediction.mat = matrix(NA, nrow=nrow(one_batch), ncol=length(q.fit))
  colnames(prediction.mat) = names(q.fit)
  for(one.mouse.idx in 1:nrow(one_batch)){
    restricted.action.idx = which(all.action.mat[,1:pulse] == action_mat[one.mouse.idx,1:pulse])
    restricted.action.mat = all.action.mat[restricted.action.idx,]
    character.actions = get_labels(restricted.action.mat)
    for(one.group.idx in 1:nrow(restricted.action.mat)){
      one.mod = q.fit[[character.actions[one.group.idx]]]
      mu.vec = one.mod$coefficients
      sig.mat = vcov(one.mod)
      naidx = which(is.na(one_batch[one.mouse.idx,]))
      non.naidx = which(!is.na(one_batch[one.mouse.idx,]))
      po = length(non.naidx)
      pm = length(naidx)
      mu.o = one_batch[one.mouse.idx,non.naidx]
      mut.o = one.mod$coefficients[non.naidx]
      mu.m = one.mod$coefficients[naidx]
      sig.oo = sig.mat[1:po,1:po]
      sig.mm = sig.mat[(po+1):(pm+po), (po+1):(pm+po)]
      sig.mo = sig.mat[(po+1):(pm+po), 1:po]
      mu.m.cond = mu.m + sig.mo%*%solve(sig.oo)%*%(mu.o - mut.o)
      sig.m.cond = sig.mm - sig.mo%*%solve(sig.oo)%*%t(sig.mo)
      if(add_variability){one.prediction.vec= MASS::mvrnorm(1, mu=mu.m.cond, Sigma=sig.m.cond)}else{
        one.prediction.vec = mu.m.cond
      }
      prediction.col = which(colnames(prediction.mat)==character.actions[one.group.idx])
      prediction.mat[one.mouse.idx, prediction.col] = one.prediction.vec[length(mu.m.cond)]
    }
  }
  
  prediction.mat
}