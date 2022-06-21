library(neuralnet)
library(reshape2)
update_q = function(q.nn, all.pairs, maxtime, inc_days){
  non.na.idx= which(!is.na(all.pairs[,maxtime,'ltv']))
  if(!is.null(q.nn)){
    start.weights=q.nn$weights
  }else{start.weights=NULL}
  all.pairs.wide = as.data.frame(cbind(cbind(all.pairs[,1:inc_days,'ltv'], all.pairs[,1:inc_days,'d']), all.pairs[,1:inc_days,'p']))
  all.pairs.wide$action = apply(all.pairs[,,'d'], c(1),get.chg.idx)-15
  q.nn = neuralnet(all.pairs[non.na.idx,maxtime, 'ltv']~.,stepmax = 10000,rep=3, thresh=1,data = all.pairs.wide[non.na.idx,], startweights=start.weights,lifesign="full", hidden=c(6,4))
  q.nn
}

get.optim.plan = function(inc.pair,q.nn, all.action.mat){
  inc.pair.wide = as.data.frame(matrix(c(inc.pair[1,1:inc_days,'ltv'], inc.pair[1,1:inc_days,'d'], inc.pair[1,1:inc_days,'p']), nrow=1))
  predictions=numeric(nrow(all.action.mat))
  for(action.idx in 1:nrow(all.action.mat)){
    inc.pair.wide$action=all.action.mat[action.idx,]
    predictions[action.idx] = predict(q.nn, inc.pair.wide)
  }
  best.plan.idx= which.min(predictions)
  all.action.mat[best.plan.idx,]
}

