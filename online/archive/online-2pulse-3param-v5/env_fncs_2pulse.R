


sequence_to_state=function(one.sequence, action.vec, done, num_free_pulses, waitime_vec=c(7,7)){
  num_zeros = 60-length(one.sequence)
  num_action_zeros = num_free_pulses-length(action.vec)
  if(num_action_zeros>0){
    action.vec=c(rep(0, num_action_zeros),action.vec)
  }
  out = c(rep(0, num_zeros),one.sequence,action.vec)
  timevec=1:length(one.sequence)
  dose_vec = rep(0, length(one.sequence))
  dose_vec[15:length(dose_vec)]=1
  rtday1 = 15+waitime_vec[1]+action.vec[1]
  if(length(dose_vec)>rtday1 & (action.vec[1]>0)){
    dose_vec[rtday1:length(dose_vec)]=2
    pd1_vec= generate_pd1_stacked(rtday1, totaltime = length(dose_vec)+7)[1:length(dose_vec)]
  }else{
    pd1_vec = generate_pd1_stacked(c(), totaltime=length(dose_vec)+7)[1:length(dose_vec)]
  }
  
  fit=lm(log(one.sequence)~timevec*dose_vec+pd1_vec)
  fit2=lm(timevec~log(one.sequence)*dose_vec+pd1_vec:dose_vec)
  ss=summary(fit)
  ss2=summary(fit2)
  out.seq=log(tail(one.sequence,20))
  out = c(as.numeric(ss$coefficients[,1:2]),ss$r.squared,ss2$coefficients[,1:2],ss2$r.squared)
  out
}

get_max = function(q.fit, one.state, potential_actions){
  
  potential.input = as.data.frame(matrix(rep(one.state, length(potential_actions)), nrow=length(potential_actions), byrow=T))
  potential.input$actions = potential_actions
  potential.input$actions2 = potential_actions^2
  #potential.input = predict(q.fit[['pcaobj']], potential.input)
  #predictions=predict(q.fit[['nnet']], potential.input[,1:6])
  predictions=predict(q.fit, potential.input)
  best_action = max(predictions)
  best_action
}

get_action = function(q.fit,one.state, potential_actions){
  
  potential.input = as.data.frame(matrix(rep(one.state, length(potential_actions)), nrow=length(potential_actions), byrow=T))
  potential.input$actions = potential_actions
  potential.input$actions2=potential_actions^2
  potential.input$actions3=potential_actions^3
  #potential.input = predict(q.fit[['pcaobj']], potential.input)
  predictions=predict(q.fit, potential.input)
  best_action = which.max(predictions)
  best_action
}

replay=function(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini,nnet_size=5, gam=.99, potential_actions=1:10){
  minisize=nrow(state_mat_mini)
  targets=rep(NA, minisize)
  for(idx in 1:minisize){
    target = rewards_mini[idx]
   
    targets[idx] = target
    
  }
  inputs=as.data.frame(state_mat_mini)
  inputs$actions=actions_mini
  inputs$actions2=actions_mini^2
  inputs$targets=targets
  
  q.fit = neuralnet(formula = targets~(.) , data=inputs, hidden = c(30), threshold = 100000,
                    stepmax = 1, rep = 1, startweights = q.fit$weights,
                    learningrate.limit = NULL, learningrate.factor = list(minus = 0.5,
                                                                          plus = 1.2), learningrate = NULL, lifesign = "none",
                    lifesign.step = 1000, algorithm = "rprop+", err.fct = "sse",
                    act.fct = "logistic", linear.output = TRUE, exclude = NULL,
                    constant.weights = NULL, likelihood = FALSE)

  #q.fit = nnet::nnet(targets~., data=inputs, size=30, Wts = q.fit$wts, maxit=1)
  q.fit

}