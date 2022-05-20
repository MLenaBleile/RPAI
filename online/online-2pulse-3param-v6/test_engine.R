setwd("~/RPAI/online/online-2pulse-3param-v6")
source("data_generation_fncs.R")
source("env_fncs_2pulse.R")

#source("online_engine.R")
set.seed(1998)
max_epochs=1000
parameter_mat = make_parameter_mat(max_epochs+500)

total_time=40
total_days=total_time
num_pulses=2
num_free_pulses=num_pulses-1
state_size = 31
bellmann_error = rep(NA, max_epochs)
waitime_vec=rep(9, num_pulses)
state_mat = matrix(NA, nrow = max_epochs, ncol=state_size)
nextstate_mat = matrix(NA, nrow=max_epochs, ncol=state_size)
potential_actions = 1:5
action_size=length(potential_actions)
actions = rep(NA, max_epochs)
dones = rep(NA, max_epochs)
eps=1
eps.vec=c(eps)
eps_decay=1
epsilon_min=.001
minibatch_size=100
epoch=1
source("online_engine.R")
#q.fit=readRDS("model.rds")

test_indices=1:50
num_mice = length(test_indices)
optimal_actions=rep(NA, num_mice)
data_mat = matrix(NA, ncol=20, nrow=num_mice)
true.mins = rep(NA, num_mice)
agent.outcomes=rep(NA, num_mice)
reference_days=c(1,potential_actions+waitime_vec[1],20)
references = c(paste("day", reference_days), "random")
reference.outcomes=matrix(NA, nrow=num_mice, ncol=length(references))
colnames(reference.outcomes) = references
selected_actions=c()
current.time=15-2+waitime_vec[1]

for(mouse in 1:num_mice){
  parameter_vec= parameter_mat[mouse,]
  one.reference = generate_one_counterfactualset(parameter_mat[mouse,],total_days = total_days+5, wait_time = 0, potential_actions = 1:20)
  one.optima= which.min(one.reference[,total_days])
  optimal_actions[mouse]= as.numeric(one.optima)
  true.mins[mouse]= log(min(one.reference[,total_days]))
  for(refday.idx in 1:length(reference_days)){
    reference_day = reference_days[refday.idx]
    #day15.outcomes[mouse] = one.reference[reference_day-6,total_days]
    reference.outcomes[mouse, refday.idx]=log(one.reference[reference_day,total_days])
  }
  one.random.action =sample(potential_actions+waitime_vec[1],1)
  reference.outcomes[mouse,'random'] = log(one.reference[one.random.action, total_days])
  parameter_vec=parameter_mat[mouse+501,]
  one.sequence = generate_one(radiation_days=c(),parameter_vec=parameter_vec,maxtime = current.time)
  one.state=sequence_to_state(one.sequence = as.numeric(one.sequence), num_free_pulses = num_free_pulses, action.vec=c())
  one.action = get_action(q.fit,one.state, potential_actions = potential_actions)+waitime_vec[1]
  selected_actions=c(selected_actions, one.action)
  one.agent.treated.sequence = generate_one(one.action+15, parameter_vec=parameter_vec, maxtime=total_days)
  agent.outcomes[mouse] = log(one.agent.treated.sequence[total_days])
}




colors = rainbow(length(references))
deltaref = reference.outcomes[,1]
plot(density(deltaref-agent.outcomes), xlab="final ltv reference - final ltv agent",type="n",xlim=c(-.2,.5), main="2 pulse performance", col="red", ylim=c(0,25))
for(refday.idx in 1:length(references)){
  deltaref = reference.outcomes[,refday.idx]
  lines(density(deltaref - agent.outcomes), col=colors[refday.idx])
  cohen = mean(deltaref-agent.outcomes)
  cat("mean: ", cohen,"\n")
  #cat("minimum:",min(deltaref-deltamodel),"\n")
  #cat("maximum:",max(deltaref-deltamodel),"\n")
  effect.sizes[refday.idx] = cohen
  #cat(deltamodel[deltamodel>deltaref])
}
legend("topright", legend=c(paste("day", reference_days), "random"), pch=16, col=c(colors, "purple"))
abline(v=0, lty=2)

print(table(selected_actions, optimal_actions))
