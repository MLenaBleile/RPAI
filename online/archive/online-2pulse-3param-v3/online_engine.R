setwd("~/Documents/Dissertation/RPAI/online/online-2pulse-3param-v3")
source("data_generation_fncs.R")
source("env_fncs_2pulse.R")

set.seed(1998)
max_epochs=8000
parameter_mat = make_parameter_mat(max_epochs+100)

total_time=40
total_days=total_time
num_pulses=2
num_free_pulses=num_pulses-1
state_size = 31
bellmann_error = rep(NA, max_epochs)
waitime_vec=rep(7, num_pulses)
state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
potential_actions = 1:10
action_size=length(potential_actions)
actions = rep(NA, max_epochs)
dones = rep(NA, max_epochs)
eps=1
eps.vec=c(eps)
eps_decay=.9995
epsilon_min=.001
minibatch_size=1000
epoch=1



rewards = c()
mouse=0
epoch=1

random.init = as.data.frame(matrix(rnorm(length(state_mat)), nrow=nrow(state_mat), ncol=state_size))
random.init$actions=sample(potential_actions,nrow(state_mat), replace=T)
random.init$actions2 = random.init$actions^2

#q.fit = lm(rnorm(nrow(state_mat))~(.)+(.)*actions+actions:actions, data=random.init)

q.fit = caret::pcaNNet(rnorm(nrow(state_mat))~(.), data=random.init, size=10)
israndom=rep(T, max_epochs)
while(epoch < max_epochs){
  mouse=mouse+1
  action.vec=c()
  current.time=15-2+waitime_vec[1]
  parameter_vec=parameter_mat[mouse+101,]
  one.sequence = generate_one(radiation_days=action.vec,parameter_vec=parameter_vec,maxtime = current.time)
  for(pulse in 1:(num_free_pulses)){
    one.state = sequence_to_state(one.sequence, action.vec, done, num_free_pulses=num_free_pulses)
    done = (pulse==num_free_pulses)
    dones[epoch]=done
    state_mat[epoch,]=one.state
    if(runif(1)>eps&(epoch>(minibatch_size+1))){one.action = get_action(q.fit, one.state, potential_actions = potential_actions); israndom[epoch]=F}else{one.action=sample(potential_actions,1)}
    action.vec=c(action.vec, one.action)
    actions[epoch]=one.action
    radiation_days = cumsum(action.vec+waitime_vec[1:length(action.vec)])+15
    
    if(done){
        current.time = total_time
        }else{
      current.time = current.time+one.action+waitime_vec[pulse+1]
    }
    one.next.sequence = generate_one(radiation_days, parameter_vec = parameter_vec, maxtime=current.time)
    one.next.state = sequence_to_state(one.sequence,action.vec,done, num_free_pulses=num_free_pulses)
    one.reward = -log(one.next.sequence[length(one.next.sequence)])/7
    rewards = c(rewards, one.reward)
    nextstate_mat[epoch,] = one.next.state
    
      if(epoch>=minibatch_size & (epoch%%minibatch_size==0)){
        minibatch_idx = 1:epoch
        state_mat_mini = state_mat[minibatch_idx,]
        nextstate_mat_mini=nextstate_mat[minibatch_idx,]
        actions_mini = actions[minibatch_idx]
        rewards_mini = rewards[minibatch_idx]
        
        dones_mini = dones[minibatch_idx]
        rewards_mini[dones_mini==0] = rewards_mini[dones_mini==1]
        
        q.fit = replay(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini)
        
        
      }
    active.eps.decay = eps_decay*(epoch>=minibatch_size) + 1*(epoch<minibatch_size)
    eps = eps*active.eps.decay
    eps.vec=c(eps.vec, eps)
    epoch = epoch+1
    cat("\n epoch:", epoch, "reward:",one.reward)
    one.sequence=one.next.sequence
  }
  # if(eps<=epsilon_min){
  #   cat("\n finished due to small epsilon")
  #   break
  # }
}

dones[epoch]=1
inputs=as.data.frame(state_mat)
inputs$actions=actions
inputs$actions2=actions^2
inputs$actions3=actions^3
train_idx = setdiff(1:4999,1:100)
q.fit = caret::pcaNNet(rewards[train_idx]~(.)^2, data=inputs[train_idx,], size=30,linout=T, scale=T, maxit=50000)
