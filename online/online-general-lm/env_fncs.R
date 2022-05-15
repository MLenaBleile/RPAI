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
    one_batch[idx,,] = generate_one(cumsum(as.numeric(action_mat[idx,]))+15, parameter_mat[idx,], current.time)
  }
  one_batch
}

get_labels = function(all.action.mat){
  char.mat = as.matrix(apply(all.action.mat,c(2),as.character))
  char.labs = char.mat[,1]
  if(ncol(char.mat)>1){
  for(jj in 2:ncol(char.mat)){
    char.labs = paste(char.labs, char.mat[,jj], sep=".")
  }}
  char.labs
}


update_model = function(all_data, params){
  oo =optim(params, obj.fnc, all_data=all_data, method="BFGS")
  oo$par
}

get_predictions = function(params, one.pair, all.action.mat){
  starting.parms = runif(3)
  names(starting.parms) = names(params)[1:3]
  one.optim = optim(starting.parms, get_loss_inc, one.pair=one.pair, method="BFGS")
  indiv.par = one.optim$par
  indiv.par[['beta.d']]=params[['beta.d']]
  indiv.par[['beta.d2']]=params[['beta.d2']]
  print(indiv.par)
  #indiv.par[['rho3']]=params[['rho3']]
  prediction.vec=rep(NA, nrow(all.action.mat))

    for(one.group.idx in 1:nrow(all.action.mat)){
      d2 = all.action.mat[one.group.idx,]
      ltv.seq = generate_one_lm(params=indiv.par, d2=d2+15)
      prediction.vec[one.group.idx] = ltv.seq[length(ltv.seq)]
    }

  prediction.vec
}

get.optim.plan = function(parameter_vec,maxtime, all.action.mat){
  ftv.vec= rep(NA, nrow(all.action.mat))
  for(j in 1:nrow(all.action.mat)){
    one.action.vec = all.action.mat[j,]
    one.sequence = generate_one(radiation_days=c(cumsum(one.action.vec)+15), parameter_vec=parameter_vec, maxtime=maxtime+10)[1,,'ltv']
    ftv.vec[j] = one.sequence[maxtime]
  }
  best.action.vec=all.action.mat[which.min(ftv.vec),]
  best.action.vec
}