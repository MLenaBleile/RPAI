library(readr)
all.round1 <- read_csv("~/Data/Round1_log.csv")

one.group = all.round1[all.round1$group == "10Gyd0d10",]
one.group.pairs = array(dim=c(nrow(one.group), 41,4))
dimnames(one.group.pairs) = list(NULL, NULL, c("ltv","d","p",'day'))
obs.days =c(10,13,17,20,23,26,29,32,37,39,41)
obs.days2=obs.days^2
obs.days3=obs.days^3
rteffects = array(0, dim=c(total_mice, length(test_effects)+1, num_free_pulses))
dimnames(rteffects) = list(1:total_mice,c(names(test_effects), 'action'), as.character(1:num_free_pulses+1))
names(dimnames(rteffects)) = c("animal", "effect", "pulse")
prediction.vec = numeric(nrow(one.group))
eval.time=32
for(mouse.idx in 1:nrow(one.group)){
  one.group.pairs[mouse.idx,,'day'] = 1:41
  dose.vec = rep(0,41)
  dose.vec[c(15,25)]=1
  pd1.vec=generate_pd1_stacked(c(25),41)
  one.group.pairs[mouse.idx,,'d'] = cumsum(dose.vec)
  one.group.pairs[mouse.idx,,'p'] = pd1.vec[1:41]
  one.real.ltv=as.numeric(one.group[mouse.idx,-c(1,2,3,4)])
  one.interpolate.fit = lm(one.real.ltv~obs.days+ obs.days2+obs.days3)
  one.group.pairs[mouse.idx,,'ltv'] = predict(one.interpolate.fit,data.frame(obs.days=1:41, obs.days2=(1:41)^2, obs.days3=(1:41)^3))
  one.group.pairs[mouse.idx,obs.days,'ltv'] = one.real.ltv
  one.pair= one.group.pairs[mouse.idx,,]
  one.effectset = get_rteffects(one.pair, num_free_pulses = 1, total_time = inc.time)
  rteffects[one.mouse,1:length(test_effects),pp] = one.effectset
  one.in.vec = one.effectset
  one.in = matrix(one.in.vec, byrow=T, nrow=1)
  colnames(one.in) = names(one.in.vec)
  one.in.df = as.data.frame(t(predict(pca.obj, as.data.frame(one.in))[,1:num_pc]))
  one.in.df$act=as.character(10)
  prediction=predict(fit2, one.in.df)
  prediction.vec[mouse.idx] = prediction 
}

print(prediction.vec)

