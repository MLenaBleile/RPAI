###


sum.exp.mean = function(d0,d1, maxtime, rho1, prop.factor, beta0, beta1, gg){
  rho2=prop.factor*rho1
  first.decay.term = rho1*loglinear.growth(d1, beta0, beta1)*tumor.decay(maxtime, radiation_days[1], gg)
  if(maxtime>d1){
    
    fct.mean.exp.1 = first.decay.term + (1-rho1)*loglinear.growth(d1, beta0, beta1)
    fct.mean.exp = rho2*fct.mean.exp.1*tumor.decay(maxtime,d1,gg)+(1-rho2)*(first.decay.term + (1-rho1)*loglinear.growth(maxtime, beta0, beta1))
    
  }else{
    fct.mean.exp = first.decay.term + (1-rho1)*loglinear.growth(maxtime, beta0, beta1)
  }
  log(fct.mean.exp)
}
#sum.exp.obj.fnc

optimize.predictive.model = function(all.pairs, selected_actions){
  df.molten=reshape2::melt(all.pairs, names(dimnames(all.pairs)))
  df.long = reshape2::dcast(df.molten,animal+time~effect, value.var="value")
  df.long$d0=15
  df.long$animal = as.integer(df.long$animal)
  df.long$d1 = selected_actions[df.long$animal]
  df.long$maxtime=maxtime
  attach(df.long)
  ####fit sexp model here
  modnlme1 <- nlme::nlme(ltv ~ sum.exp.mean(d0,d1, maxtime, rho1, prop.factor, beta0, beta1, gg),fixed=c(rho1~1, prop.factor~1, beta0~1, beta1~1, gg~1), data = df.long, start=rep(1,5))
  
  estimated.par
}

get.predicted.optim.plan = function(predictive.params,maxtime, all.action.mat){
  all.cf.vec = predict.sexp(predictive.params, maxtime, all.action.mat)
  min.idx=which.min(all.cf.vec)
  one.action= all.action.mat[min.idx,]
  one.action
}

predict.sexp = function(predictive.params, maxtime, one.action){
  p1.vec = generate_pd1_stacked(all.action.mat$Var1+15, maxtime)
  dose_vec = rep(0,maxtime)
  dose_vec[one.action+15] = 10
  dose_vec[15] = 10
  rho.inv.logit = predictive.params['alphaio']*p1.vec + predictive.params['alphart']*dose_vec
  
}