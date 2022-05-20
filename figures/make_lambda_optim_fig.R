nn=500
lambda.vec=seq(0.01,0.5, length.out=nn)
names(lambda.vec) = rep("lambda",nn)
optims=numeric(nn)
all.action.mat=make_potential_action_mat(5:15,1)
for(jj in 1:nn){
  optims[jj] = get.optim.plan(parameter_vec = c(lambda.vec[jj]), maxtime=40, all.action.mat = all.action.mat
                              )
}
par(mfrow=c(1,1))
beanplot::beanplot(lambda.vec~optims, xlab="optimal day spacing", ylab=expression(lambda))
