setwd("~/RPAI/online/online-4pulse-3param-v2")
source("data_generation_fncs.R")
source("env_fncs.R")

set.seed(1998)
max_epochs=3200
test_num=500
parameter_mat = make_parameter_mat(max_epochs+test_num)

total_time=120
num_pulses=4
num_free_pulses=num_pulses-1
state_size = total_time+11+num_pulses
#state_size = total_time+num_pulses
#state_size=32

waitime_vec=rep(8, num_free_pulses)
state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
potential_actions = 1:13
action_size=length(potential_actions)


bellmann_error = rep(NA, max_epochs)
actions = rep(NA, max_epochs)
dones = rep(NA, max_epochs)
eps=1
eps.vec=c(eps)
eps_decay=.999
burn_in=100
epsilon_min=.001
minibatch_size=100
epoch=1




rewards = c()
mouse=1
epoch=1
mouse2=0

random.init = as.data.frame(matrix(rnorm(length(state_mat)), nrow=nrow(state_mat), ncol=state_size))
random.init$actions=sample(potential_actions,nrow(state_mat), replace=T)
random.init$actions2 = random.init$actions^2

#q.fit = lm(rnorm(nrow(state_mat))~(.)+(.)*actions+actions:actions, data=random.init)
library(neuralnet)
q.fit = neuralnet(formula = V1~(.) , data=random.init, hidden = c(10,6,4), threshold = 10000,stepmax = 2, rep = 2)
israndom=rep(T, max_epochs)
while(epoch < max_epochs){

  action.vec=c()
  current.time=15-2+waitime_vec[1]
  parameter_vec=parameter_mat[mouse+test_num+1,]
  one.sequence = generate_one(radiation_days=action.vec,parameter_vec=parameter_vec,maxtime = current.time)
  for(pulse in 1:(num_free_pulses)){
    one.state = sequence_to_state(one.sequence, action.vec, done, num_free_pulses=num_free_pulses, total_time=total_time)
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
    one.reward = -(log(one.next.sequence[length(one.next.sequence)]) +16)/6
    rewards = c(rewards, one.reward)
    # if(done){
    #   rewards[1:(epoch) & !(dones[1:epoch])] = rewards[1:(epoch) & dones[1:epoch]]
    # }
    nextstate_mat[epoch,] = one.next.state
    
      if(epoch>=burn_in & done){
        minibatch_idx = sample(1:epoch, minibatch_size)
        state_mat_mini = state_mat[minibatch_idx,]
        nextstate_mat_mini=nextstate_mat[minibatch_idx,]
        actions_mini = actions[minibatch_idx]
        rewards_mini = rewards[minibatch_idx]
        
        dones_mini = dones[minibatch_idx]
        #rewards_mini[dones_mini==0] = 0
        
        q.fit = replay(q.fit,state_mat_mini,actions_mini,nextstate_mat_mini, rewards_mini, dones_mini)
        
        if(eps>epsilon_min){eps = eps*eps_decay}
        
      }
    #active.eps.decay = eps_decay*(epoch>=burn_in) + 1*(epoch<burn_in)
    
    
    eps.vec=c(eps.vec, eps)
    epoch = epoch+1
    cat("\n epoch:", epoch, "reward:",one.reward)
    one.sequence=one.next.sequence
  }
}

 minibatch_idx = 1:(max_epochs-1)
 inputs= as.data.frame(state_mat[minibatch_idx,])
 inputs$actions=actions[minibatch_idx]
 inputs$actions2=actions[minibatch_idx]^2
 inputs$targett=rewards[minibatch_idx]
q.fit = neuralnet::neuralnet(formula = targett~.,stepmax=10000,thresh=.5,rep=5,data=inputs, hidden = c(10,6,4), lifesign='full')
  
