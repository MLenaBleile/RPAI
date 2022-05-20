x11()
one.sequence=generate_one(c(),parameter_vec, maxtime=maxtime)[1,,'ltv']
par(mfrow=c(2,2))
plot(one.sequence, xlab="days", ylab="ltv", type="l")
colors=rainbow(10)
labels=matrix(nrow=10, ncol=3)
for(jj in 1:10){
  one.sequence=generate_one(c(),minibatch_parameters[jj,], maxtime=maxtime)[1,,'ltv']
  lines(one.sequence, xlab="days", ylab="ltv", type="l", col=colors[jj])
  labels[jj,]=paste(paste(random_params,"="),as.character(round(minibatch_parameters[jj,random_params], digits=3)))
}

one.sequence=exp(generate_one(c(),parameter_vec, maxtime=maxtime)[1,,'ltv'])
plot(one.sequence, xlab="days", ylab="tv", type="l")
for(jj in 1:10){
  one.sequence=exp(generate_one(c(),minibatch_parameters[jj,], maxtime=maxtime)[1,,'ltv'])
  lines(one.sequence, col=colors[jj])
  labels[jj,]=paste(paste(random_params,"="),as.character(round(minibatch_parameters[jj,random_params], digits=3)))
}

plot(one.sequence, xlab="parameters", ylab="", type="n")
paste.all = function(str.vec){
  str.str = str.vec[1]
  for(j in 2:length(str.vec)){
    str.str=paste(str.str, str.vec[j])
  }
  str.str
}
legend("topleft", legend=apply(labels, c(1), paste.all), col=colors, pch=16)
