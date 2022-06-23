
sequence_to_state=function(one.sequence, action.vec, done, num_free_pulses, wait_time, total_time){
  num_zeros = total_time-length(one.sequence) - wait_time
  num_action_zeros = num_free_pulses-length(action.vec)
  if(num_action_zeros>0){
    action.vec=c(rep(0, num_action_zeros),action.vec)
    
  }
  out = c(rep(0, num_zeros),log(one.sequence),action.vec)
  timevec=1:length(one.sequence)
  dose_vec = rep(0, length(one.sequence))
  dose_vec[15:length(dose_vec)]=1
  rtday1 = 15+wait_time+action.vec[1]
  
  if(length(dose_vec)>rtday1 & (action.vec[1]>0)){
    dose_vec[rtday1:length(dose_vec)]=2
    pd1_vec= generate_pd1_stacked(rtday1, totaltime = length(dose_vec)+7)[1:length(dose_vec)]
  }else{
    pd1_vec = generate_pd1_stacked(c(), totaltime=length(dose_vec)+7)[1:length(dose_vec)]
  }
  
  fit=lm(log(one.sequence)~timevec*dose_vec*pd1_vec)
  ss=summary(fit)
  dose_rescale=invlogit((max(dose_vec)+.5)/3)
  time_rescale = invlogit((length(one.sequence)+.5)/61)
  #out = c(as.numeric(ss$coefficients[,1]),time_rescale, ss$r.squared, dose_rescale)
  out = c(out, time_rescale,dose_rescale,as.numeric(ss$coefficients[,1]))
}

get_action = function(q.fit,one.state, potential_actions){
  
  potential.input = as.data.frame(matrix(rep(one.state, length(potential_actions)), nrow=length(potential_actions), byrow=T))
  #potential.input = as.data.frame(stats::prcomp(potential.input, scale=T)$x)
  potential.input$actions = as.character(potential_actions)
  #potential.input$actions2=potential_actions^2
  #colnames(potential.input) = make.names(colnames(potential.input))
  #potential.input = predict(q.fit[['pcaobj']], potential.input)
  
  form = ~(.)*actions*(.)
  modmat = as.data.frame(model.matrix(form, data = potential.input))
  colnames(modmat) = make.names(colnames(modmat))
  predictions=predict(q.fit, modmat)
  best_action = which.max(predictions)
  best_action
}

replay=function(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini,nnet_size=5,stepmax=1, reps=1,threshold=100000, gam=.99, potential_actions=1:10){
  minisize=nrow(state_mat_mini)
  targets=rep(NA, minisize)
  for(idx in 1:minisize){
    target = rewards_mini[idx]
    #if(FALSE){target=target + gam*get_max(q.fit, one.state=nextstate_mat_mini[idx,], potential_actions = potential_actions)}
    targets[idx] = target
    
  }
  inputs=as.data.frame(state_mat_mini)
  #pca.obj = stats::prcomp(inputs)
  #inputs = as.data.frame(pca.obj$x)
  inputs$actions=as.character(actions_mini)
  #inputs$actions2=actions_mini^2
  
  
  #q.fit = nnet::nnet(targets~(.)+(.)*actions, data=scale(inputs), size=nnet_size)
  
  
  inputs$targets=targets
  form = targets~(.)*actions
  modmat = as.data.frame(model.matrix(form, data = inputs))
  modmat$targets=targets
  colnames(modmat) = make.names(colnames(modmat))

  q.fit = neuralnet(targets~(.), data=modmat, hidden = c(10,6,4), threshold = threshold,
                    stepmax = stepmax, rep = reps, startweights = q.fit$weights,
                    learningrate.limit = NULL, learningrate.factor = list(minus = 0.5,
                                                                          plus = 1.2), learningrate = NULL, lifesign = "full",
                    lifesign.step = 1000, algorithm = "rprop+", err.fct = "sse",
                    act.fct = "logistic", linear.output = TRUE, exclude = NULL,
                    constant.weights = NULL, likelihood = FALSE)
  
  q.fit
}