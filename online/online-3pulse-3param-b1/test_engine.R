

# #q.fit=readRDS("model.rds")
# test_num=50
# 
# 
# train_num = dim(all_data)[1]
# animal.idx=dim(all_data)[1]
# test_idx = (train_num+1):(train_num+test_num)
# 
# while(animal.idx<(test_num+train_num)){
#   cat("animal:", animal.idx,"\n")
#   test_parameter_mat = make_test_parameter_mat(5)
#   one.pair = generate_one(c(), parameter_mat[1,], inc.time)
#   
#   num.extra = -dim(one.pair)[2] + dim(all_data)[2]
#   one.pair.extended = abind::abind(one.pair, array(NA, dim=c(1, num.extra,4)), along=2)
#   dimnames(one.pair.extended) = list(animal.idx+1, 1:total_time,c("ltv","d","p","day"))
#   q.fit = update_model(q.fit,all_data, one.pair.extended, all_actions = all_actions, all.action.mat)
#   
#   prediction.vec = get_predictions(q.fit, one.pair.extended, all.action.mat, add_variability = F)
#   best_action = names(q.fit)[which.min(prediction.vec)]
#   action_mat = as.matrix(as.numeric(best_action), ncol=num_free_pulses)
#   parameter_mat=t(as.matrix(parameter_mat[1,], nrow=1, ncol=8))
#   one_updated_batch = generate_one_batch(parameter_mat, action_mat, current.time = total_time)
#   all_data = abind::abind(all_data, one_updated_batch,along=1)
#   
#   all_parameters = rbind(all_parameters, parameter_mat)
#   all_actions = rbind(all_actions, action_mat)
#   animal.idx = nrow(all_data)
#   dimnames(all_data) = list(1:animal.idx, 1:total_time, c("ltv","d","p","day"))
# }


test_num=length(test_mice)
all_test_data = one_batch[test_mice,,]
all_test_parameters = parameter_mat[test_mice,]
all_test_actions = act.mat[test_mice]

reference_days = c(1,potential_actions,20)
optim_actions=c()
ref.outcome.mat = matrix(NA, nrow=test_num, ncol=length(reference_days)+1)
references=c(paste("day",reference_days),"random")
colnames(ref.outcome.mat)=references
test.id.sequential=1
for(test.id in test_mice ){
  
  for(refday.idx in 1:length(reference_days)){
    action_mat = as.matrix(rep(reference_days[refday.idx],num_free_pulses), ncol=num_free_pulses)
    test_parameter_mat=t(as.matrix(all_test_parameters[test.id.sequential,], nrow=1, ncol=8))
    one_ref_sequence = generate_one_batch(test_parameter_mat, t(action_mat), current.time = total_time)
    ref.outcome.mat[test.id.sequential,refday.idx] = one_ref_sequence[,total_time,'ltv']
  }
  action_mat = as.matrix(all.action.mat[sample(1:nrow(all.action.mat),1),])
  test_parameter_mat=t(as.matrix(all_test_parameters[test.id.sequential,], nrow=1, ncol=8))
  one_ref_sequence = generate_one_batch(test_parameter_mat, action_mat, current.time = total_time)
  ref.outcome.mat[test.id.sequential,'random'] = one_ref_sequence[,total_time,'ltv']
  #optim_actions=c(optim_actions, get.optim.plan(test_parameter_mat[1,], maxtime=total_time, all.action.mat = as.matrix(all.action.mat)))
  test.id.sequential= test.id.sequential+1
  }


adj=1
loss.mat = matrix(NA, nrow=test_num, ncol=ncol(ref.outcome.mat))
plot(density(ref.outcome.mat[,"random"]-all_test_data[, total_time,'ltv'],adjust=adj),xlim=c(-.2,.5),main="3 pulse performance", type="n", xlab="ltv reference-ltv agent")
colourss= rainbow(ncol(ref.outcome.mat))
loss.means = rep(NA, length(references))
loss.medians = rep(NA, length(references))
loss.mins = rep(NA, length(references))
loss.maxes=rep(NA, length(references))
for(refday.idx in 1:ncol(ref.outcome.mat)){
  loss=ref.outcome.mat[,refday.idx]-all_test_data[, total_time,'ltv']
  lines(density(loss, adjust = adj), col=colourss[refday.idx])
  loss.mat[,refday.idx]=loss
  loss.means[refday.idx] = mean(loss)
  loss.medians[refday.idx] = median(loss)
  loss.mins[refday.idx] = min(loss)
  loss.maxes[refday.idx] = max(loss)
}
legend("topright", legend=paste(colnames(ref.outcome.mat), round(loss.means,digits=4), sep=": "), col=colourss, pch=16)
abline(v=0, lty=2)
#print(table(all_test_actions, optim_actions))

colnames(loss.mat) = colnames(ref.outcome.mat)
loss.df = stack(as.data.frame(loss.mat))
result = data.frame(means=loss.means, medians=loss.medians, mins=loss.mins, maxes=loss.maxes)
row.names(result) = colnames(ref.outcome.mat)
print(result)
