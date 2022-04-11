make_potential_action_mat = function(potential_actions, num_free_pulses){
  expand.grid(rep(list(potential_actions), num_free_pulses))
}


generate_one_initial_batch = function(parameter_mat,inc.time, total_time){
  minibatch_size = nrow(parameter_mat)
  one_batch = array(NA, dim=c(minibatch_size, total_time, 4))
  dimnames(one_batch) = list(NULL, NULL, c("ltv","d","p","day"))
  for(idx in 1:minibatch_size){
    one_batch[idx,1:inc.time,] = generate_one(c(), parameter_mat[idx,], inc.time)

  }
  one_batch
}

generate_one_batch = function(parameter_mat, action_mat, current.time){
  minibatch_size = nrow(parameter_mat)
  one_batch = array(NA, dim=c(minibatch_size, total_time, 4))
  dimnames(one_batch) = list(NULL, NULL, c("ltv","d","p","day"))
  for(idx in 1:minibatch_size){
    one_batch[idx,,] = generate_one(cumsum(action_mat[idx,])+15, parameter_mat[idx,], current.time)
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

update_model = function(q.fit,all_data, one.pair, action_mat){
  character.actions = get_labels(action_mat)
  unique.trt = unique(character.actions)
  for(one.group.idx in 1:length(unique.trt)){
  one.group.mice = which(character.actions==unique.trt[one.group.idx])
  num.extra = -dim(one.pair)[2] + dim(all_data)[2]
  one.pair.extended = abind::abind(one.pair, array(NA, dim=c(1, num.extra,4)), along=2)
  one_ds = abind::abind(all_data[one.group.mice,,],one.pair.extended, along=1)

  one_ds.long = plyr::adply(one_ds, c(1))

  one.model = glm(ltv~day +day:d+ p:d, data=one_ds.long)
  q.fit[[unique.trt[one.group.idx]]] = one.model
  }
  q.fit
}

get_predictions = function(q.fit, one.pair.extended, all.action.mat, add_variability=T){
  character.actions = get_labels(all.action.mat)
  prediction.vec = rep(NA, length(character.actions))

    for(one.group.idx in 1:nrow(all.action.mat)){
      one.mod = q.fit[[character.actions[one.group.idx]]]
      mu.vec = one.mod$coefficients
      sig.mat = vcov(one.mod)
      naidx = which(is.na(one.pair.extended[,,'ltv']))
      non.naidx = which(is.na(one.pair.extended[,,'ltv']))
      po = length(non.naidx)
      pm = length(naidx)
      mu = one.mod$coefficients
      mu[-c(1)] = mu[-c(1)]+mu[1]
      mu.o = one_batch[one.mouse.idx,non.naidx]
      mut.o = mu[non.naidx]
      mu.m = mu[naidx]
      sig.oo = sig.mat[1:po,1:po]
      sig.mm = sig.mat[(po+1):(pm+po), (po+1):(pm+po)]
      sig.mo = sig.mat[(po+1):(pm+po), 1:po]
      mu.m.cond = mu.m + sig.mo%*%solve(sig.oo)%*%(mu.o - mut.o)
      sig.m.cond = sig.mm - sig.mo%*%solve(sig.oo)%*%t(sig.mo)
      if(add_variability){one.prediction.vec= MASS::mvrnorm(1, mu=mu.m.cond, Sigma=sig.m.cond)}else{
        one.prediction.vec = mu.m.cond
      }
      prediction.mat[one.mouse.idx, one.group.idx] = mean(one.prediction.vec[(length(mu.m.cond)-3):length(mu.m.cond)])
    }
  
  colnames(prediction.mat) = character.actions
  prediction.mat
}