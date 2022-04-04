setwd("~/RPAI/online/online-5pulse-working")
source("data_generation_fncs.R")
source("env_fncs.R")
library(neuralnet)

set.seed(1998)
max_epochs=5000
test_num=500
parameter_mat = make_parameter_mat(max_epochs+test_num)

total_time=120
num_pulses=5
num_free_pulses=num_pulses-1
state_size = total_time+11+num_pulses
#state_size = total_time+num_pulses
#state_size=12

waitime_vec=rep(9, num_free_pulses)
state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
potential_actions = 1:15
action_size=length(potential_actions)


bellmann_error = rep(NA, max_epochs)
actions = rep(NA, max_epochs)
dones = rep(NA, max_epochs)
eps=1
eps.vec=c(eps)
eps_decay=.998
burn_in=150
epsilon_min=.001
minibatch_size=150
epoch=1


rewards = c()
mouse=1
epoch=1
mouse2=0
refit=T 

random.init = as.data.frame(matrix(rnorm(length(state_mat)), nrow=nrow(state_mat), ncol=state_size))
random.init$actions=sample(potential_actions,nrow(state_mat), replace=T)
random.init$actions2 = random.init$actions^2
random.init$targets = rnorm(nrow(random.init))

#q.fit = lm(rnorm(nrow(state_mat))~(.)+(.)*actions+actions:actions, data=random.init)
#library(neuralnet)
q.fit = neuralnet::neuralnet(formula = targets~(.) , data=random.init, hidden = c(10,6), threshold = 10000,stepmax = 2, rep = 2)
israndom=rep(T, max_epochs)

while(epoch < max_epochs){

  action.vec=c()
  current.time=15-2+waitime_vec[1]
  parameter_vec=parameter_mat[mouse+1+test_num,]
  one.sequence = generate_one(radiation_days=action.vec,parameter_vec=parameter_vec,maxtime = current.time)
  for(pulse in 1:(num_free_pulses)){
    one.state = sequence_to_state(one.sequence, action.vec, done, num_free_pulses=num_free_pulses, waitime_vec = waitime_vec, total_time=total_time)
    done = (pulse==num_free_pulses)
    dones[epoch]=done
    state_mat[epoch,]=one.state
    if(runif(1)>eps&(epoch>(minibatch_size+1))){one.action = get_action(q.fit, one.state, potential_actions = potential_actions); israndom[epoch]=F}else{one.action=sample(potential_actions,1)}
    action.vec=c(action.vec, one.action)
    actions[epoch]=one.action
    radiation_days = cumsum(action.vec+waitime_vec[1:length(action.vec)])+15
    
    if(done){
        current.time = total_time
        mouse=mouse+1
        }else{
      current.time = current.time+one.action+waitime_vec[pulse+1]
    }
    one.next.sequence = generate_one(radiation_days, parameter_vec = parameter_vec, maxtime=current.time)
    one.next.state = sequence_to_state(one.sequence,action.vec,done, num_free_pulses=num_free_pulses, total_time=total_time)
    one.reward = -log(one.next.sequence[length(one.next.sequence)])/9
    rewards = c(rewards, one.reward)
    nextstate_mat[epoch,] = one.next.state
    if(done){
      rewards[1:(epoch) & !(dones[1:epoch])] = rewards[1:(epoch) & dones[1:epoch]]
    }
    
      if(epoch>=burn_in & (done)){
        minibatch_idx = sample(1:epoch, minibatch_size)
        state_mat_mini = state_mat[minibatch_idx,]
        nextstate_mat_mini=nextstate_mat[minibatch_idx,]
        actions_mini = actions[minibatch_idx]
        rewards_mini = rewards[minibatch_idx]
        
        dones_mini = dones[minibatch_idx]
        #rewards_mini[dones_mini==0] = rep(rewards_mini[dones_mini==1], each=num_free_pulses-1)
        
        q.fit = replay(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini, reps=1)
        
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


q.fit = replay(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini, reps=3,stepmax=1500,threshold = .1)

