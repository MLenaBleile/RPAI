
real.dat <- read.csv("~/Data/Round1_log.csv", row.names=1)
interest.groups = c("10Gyd0d10_PDL")
non.na.idx=  which(!is.na(real.dat$Day39))
group.idx = intersect(which(real.dat$group%in% interest.groups), non.na.idx)
result = real.dat[group.idx,c("group","Day39")]
days = c(10,13,17,20,23,26,29,32,37,39,41)
pred.rewards=numeric()
all.real.input = matrix(NA, length(group.idx), 23)
colnames(all.real.input) = c(paste("V",1:21, sep=""), "actions","actions2")
rownames(all.real.input) = group.idx
for(real.mouse in group.idx){
  one.real.mouse = real.dat[real.mouse,-c(1:3, which(is.na(real.dat[real.mouse,])))]
  one.days = days[1:length(one.real.mouse)]
  one.days2=one.days^2
  one.days3=one.days^3
  one.days4=one.days^4
  prediction.data = data.frame(ltv=as.numeric(one.real.mouse), days=one.days, days2=one.days2, days3=one.days3, days4=one.days4)
  one.fit = lm(ltv~., data=prediction.data)
  one.real.sequence.1 = predict(one.fit, data.frame(days=1:20, days2=(1:20)^2,days3=c(1:20)^3, days4=(1:20)^4))
  one.real.state = sequence_to_state(c(1,one.real.sequence.1[1:19]), action.vec=c(3), done=T, num_free_pulses = 1,eval.time=39)
  one.real.state=c(one.real.state, 3, 3^2)
  names(one.real.state) = c(paste("V",1:21, sep=""), "actions","actions2")
  all.real.input[as.character(real.mouse),] = one.real.state
  one.pred.df = as.data.frame(t(one.real.state))
  pred.rewards[real.mouse] = predict(q.fit, one.pred.df)
}
pred.rewards=pred.rewards[group.idx]
pred.rewards -  real.dat$Day39[group.idx]

q.fit.calibrate = neuralnet(result$Day39~(.), data=all.real.input, hidden=c(6,4), startweights = q.fit$weights)
