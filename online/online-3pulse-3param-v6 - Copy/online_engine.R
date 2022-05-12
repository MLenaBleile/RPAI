setwd("~/RPAI/online/online-3pulse-3param-v6 - Copy")
  source("data_generation_fncs.R")
  source("env_fncs.R")
  
  set.seed(1998)
  max_epochs=1500
  test_num=50
  num_pulses=3
  num_free_pulses=num_pulses-1
  library(LaplacesDemon)
  parameter_mat = make_parameter_mat(max_epochs/num_free_pulses+test_num)
  
  total_time=60
  
  #state_size = total_time+11+num_pulses+1
  #state_size = total_time+num_pulses+1
  state_size=63
  #state_size=
  
  wait_time=9
  state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
  nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
  potential_actions = 1:5
  action_size=length(potential_actions)
  
  
  bellmann_error = rep(NA, max_epochs)
  actions = rep(NA, max_epochs)
  dones = rep(NA, max_epochs)
  eps=1
  eps.vec=c(eps)
  eps_decay=.995
  burn_in=100
  epsilon_min=.001
  minibatch_size=64
  epoch=1
  
  inc.days=15-2+wait_time
  
  rewards = c()
  mouse=0
  epoch=1
  used_mice =c()
  random.init = as.data.frame(matrix(rnorm(length(state_mat)), nrow=nrow(state_mat), ncol=state_size))
  
  random.init$actions=as.character(sample(potential_actions,nrow(state_mat), replace=T))
  #random.init$actions2 = random.init$actions^2
  random.init$targets = rnorm(nrow(random.init))
  #colnames(random.init) = make.names(colnames(random.init))
  
  form = targets~(.)*actions
  library(neuralnet)
  modmat = as.data.frame(model.matrix(form, data = random.init))
  colnames(modmat) = make.names(colnames(modmat))
  modmat$targets=random.init$targets
  q.fit = neuralnet(targets~(.), data=modmat, hidden=c(10,6,4), stepmax=1,rep=3, threshold = 1000000)

  
  #q.fit = neuralnet::neuralnet(formula = targets~(.) , data=random.init, hidden = c(10,6), threshold = 10000,stepmax = 2, rep = 2)
  israndom=rep(T, max_epochs)
  while(epoch < max_epochs){
    mouse=mouse+1
    action.vec=c()
    current.time=inc.days
    parameter_vec=parameter_mat[mouse,]
    used_mice = c(used_mice, mouse)
    one.sequence = generate_one(radiation_days=action.vec,parameter_vec=parameter_vec,maxtime = current.time)
    for(pulse in 1:(num_free_pulses)){
      one.state = sequence_to_state(one.sequence, action.vec, done, num_free_pulses=num_free_pulses,wait_time=wait_time, total_time=total_time)
      done = (pulse==num_free_pulses)
      dones[epoch]=done
      state_mat[epoch,]=one.state
      if(runif(1)>eps&(epoch>(minibatch_size+1))){
        one.action = get_action(q.fit, one.state,potential_actions = potential_actions); israndom[epoch]=F
        }else{one.action=sample(potential_actions,1)}
      action.vec=c(action.vec, one.action)
      actions[epoch]=one.action
      radiation_days = cumsum(action.vec+wait_time)+15
      
      if(done){
          current.time = total_time
          
          }else{
        current.time = current.time+one.action+wait_time
      }
      one.next.sequence = generate_one(radiation_days, parameter_vec = parameter_vec, maxtime=current.time)
      one.next.state = sequence_to_state(one.sequence,action.vec,done, num_free_pulses=num_free_pulses,wait_time=wait_time, total_time=total_time)
      one.reward = -log(one.next.sequence[total_time])/8
      rewards = c(rewards, one.reward)
      if(done){
        rewards[1:(epoch) & !(dones[1:epoch])] = rewards[1:(epoch) & dones[1:epoch]]
      }
      nextstate_mat[epoch,] = one.next.state

        if(epoch>=burn_in & (done)){
          minibatch_idx = sample(1:epoch,minibatch_size)
          state_mat_mini = state_mat[minibatch_idx,]
          nextstate_mat_mini=nextstate_mat[minibatch_idx,]
          actions_mini = actions[minibatch_idx]
          rewards_mini = rewards[minibatch_idx]
          
          dones_mini = dones[minibatch_idx]
          
          
          q.fit = replay(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini,reps=3, stepmax=3, threshold=1000)
          if(eps>epsilon_min){eps = eps*eps_decay}
          
        }
      eps.vec=c(eps.vec, eps)
      epoch = epoch+1
      cat("\n epoch:", epoch, "reward:",one.reward)
      
      one.sequence=one.next.sequence
    }
  }

  minibatch_idx = 1:(epoch-1)
  state_mat_mini = state_mat[minibatch_idx,]
  nextstate_mat_mini=nextstate_mat[minibatch_idx,]
  actions_mini = actions[minibatch_idx]
  rewards_mini = rewards[minibatch_idx]
  
  dones_mini = dones[minibatch_idx]
  
  
  q.fit = replay(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini,reps=5, stepmax=15000, threshold=.5)
  
  source("test_engine_multi.R")
  