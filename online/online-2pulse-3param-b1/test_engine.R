

#q.fit=readRDS("model.rds")
test_num=2*minibatch_size*50
all.action.mat = make_potential_action_mat(potential_actions = potential_actions, num_free_pulses = num_free_pulses)
parameter_mat = make_parameter_mat(minibatch_size)
action_mat = matrix(sample(potential_actions, minibatch_size, replace=T),byrow=T, nrow=minibatch_size, ncol=num_free_pulses)
one_batch = generate_one_batch(parameter_mat = parameter_mat, action_mat, current.time = total_time)

all_test_data= matrix(NA, ncol=total_time, nrow=0)
all_test_actions = matrix(NA, ncol=num_free_pulses, nrow=0)
all_test_parameters = matrix(NA, ncol=8, nrow=0)
animal.idx=nrow(all_test_data)

while(animal.idx<test_num){
  parameter_mat = make_parameter_mat(test_num)
  one_batch = generate_one_initial_batch(parameter_mat, inc.time, total_time)
  #q.fit = update_model(q.fit,all_test_data, one_batch, action_mat)
  prediction.mat = get_predictions(q.fit, one_batch, all.action.mat, add_variability = F)
  best_actions = names(q.fit)[apply(prediction.mat, c(1), which.min)]
  action_mat = as.matrix(as.numeric(best_actions), ncol=num_free_pulses)
  one_updated_batch = generate_one_batch(parameter_mat, action_mat, current.time = total_time)
  all_test_data = rbind(all_test_data, one_updated_batch)
  all_test_parameters = rbind(all_test_parameters, parameter_mat)
  all_test_actions = rbind(all_test_actions, action_mat)
  animal.idx = nrow(all_test_data)
}

reference_days = c(1,10,14,20)
ref.outcome.mat = matrix(NA, nrow=test_num, ncol=length(reference_days)+1)
references=c(paste("day",reference_days),"random")
colnames(ref.outcome.mat)=references
for(refday.idx in 1:length(reference_days)){
  ref_action_mat = matrix(reference_days[refday.idx],nrow=test_num, ncol=num_free_pulses)
  reference_batch = generate_one_batch(all_test_parameters, ref_action_mat, total_time)
  ref.outcome.mat[,refday.idx] = reference_batch[,total_time]
}
random_action_mat = matrix(sample(potential_actions, test_num*num_free_pulses, replace=T), nrow=test_num, ncol=num_free_pulses)
random_batch = generate_one_batch(all_test_parameters, random_action_mat, total_time)
ref.outcome.mat[,"random"] = random_batch[,total_time]


plot(density(ref.outcome.mat[,"random"]-all_test_data[, total_time]),xlim=c(-.1,.5),main="2 pulse performance", type="n", xlab="ltv reference-ltv agent")
colourss= rainbow(ncol(ref.outcome.mat))
loss.means = rep(NA, length(references))
for(refday.idx in 1:ncol(ref.outcome.mat)){
  loss=ref.outcome.mat[,refday.idx]-all_test_data[, total_time]
  lines(density(loss), col=colourss[refday.idx])
  loss.means[refday.idx] = mean(loss)
}
legend("topright", legend=paste(colnames(ref.outcome.mat), round(loss.means,digits=3), sep=": "), col=colourss, pch=16)
abline(v=0, lty=2)
print(table(all_test_actions))
