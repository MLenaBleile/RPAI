


sequence_to_state=function(one.sequence, action.vec, done, num_free_pulses,mouse, waitime_vec=c(7,7)){
  
  out=matrix(c(log(one.sequence), 1:length(one.sequence)),ncol=2, byrow=F)
  colnames(out) = c("ltv", "day")
  out
}

get_max = function(q.fit, one.state, potential_actions){
  
  potential.input = as.data.frame(matrix(rep(one.state, length(potential_actions)), nrow=length(potential_actions), byrow=T))
  potential.input$actions = potential_actions
  potential.input$actions2 = potential_actions^2
  predictions=predict(q.fit, potential.input)
  best_predicted_reward = max(predictions[,1])
  best_predicted_reward
}

get_action = function(q.fit,one.state, potential_actions){
  
  ##do the conditional probability thing here
  state_stack = one.state

  for(jj in 2:length(potential_actions)){
    state_stack = cbind(state_stack, state_stack)
  }
  potential.input = as.data.frame(state_stack)
  potential.input$actions = rep(potential_actions, each= nrow(one.state))
  potential.input$actions2=potential_actions^2
  potential.input$actions3=potential_actions^3
  predictions=posterior_predict(q.fit, potential.input)
  best_action = which.max(predictions)
  best_action
}

replay=function(q.fit,state_mat_mini,actions_mini, rewards_mini, last_parameters){

  inputs=as.data.frame(state_mat_mini)
  inputs$actions=actions_mini
  inputs$actions2=actions_mini^2
  inputs$rewards = rewards_mini
  
  q.fit = brms::brm(targets~ltv*actions*as.factor(day) +ltv*actions2,data=inputs, prior = )

  q.fit

}