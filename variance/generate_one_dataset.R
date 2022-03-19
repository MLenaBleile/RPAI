







###generate one dataset given noise vec
generate_one_dataset = function(noise_vec=c(), parameter_vec=c(), old_names, one.round){
set.seed(1998)

new_para_set=make_default_para_set()
one.round.simulated = as.data.frame(matrix(NA, nrow=nrow(one.round), ncol=41))
one.round.simulated$group=NA
one.round.simulated$ID = one.round$ID
names_for_loop=unique(c(names(parameter_vec), names(noise_vec)))
new_para_vec = rep(NA, length(names_for_loop))
names(new_para_vec) = names_for_loop
# if(length(parameter_vec)>0){
#   for(parameter in names_for_loop){
#     new_para_set[parameter,'best'] = 
#     #new_para_vec[parameter] = truncnorm(1,loc =new_para_set[parameter,"best"],scale=noise_vec[parameter], lwr=new_para_set[parameter,"lower"], upr = new_para_set[parameter, "upper"])
#   }
# }
num_params = length(names_for_loop)

cat("\n new para vec", new_para_vec)

  for(mouse in unique(one.round$ID)){
    group_name = as.character(one.round$group)[one.round$ID==mouse]
      for(parameter in names_for_loop){
        #new_para_set[parameter,'best'] = parameter_vec[parameter]
        if(noise_vec[parameter]==0){
          new_para_vec[parameter] = parameter_vec[parameter]
        }else if(parameter=="xxx"){
          one.loc = parameter_vec[parameter]
          one.scale = noise_vec[parameter]
          rbeta(1,shape1 = one.loc*one.scale, shape2 = (1-one.loc)*one.scale)
          new_para_vec[parameter] = runif(1)
          print(new_para_vec)
          print(one.loc)
          print(one.scale)
          #quit()
        }else{
        new_para_vec[parameter] = truncnorm(1,loc =parameter_vec[parameter],scale=noise_vec[parameter], lwr=new_para_set[parameter,"lower"], upr = new_para_set[parameter, "upper"])
        
        }
      }
    
    
  
    is.pd1group = (length(unlist(strsplit(as.character(group_name),"_")))-1 ==1)
  onegp.split = suppressWarnings(as.numeric(unlist(strsplit(as.character(group_name),"d"))))
  rt.days = onegp.split[!is.na(onegp.split)]+15
  #rt.days = rt.days[rt.days!=15]
  one.seq = one.round[one.round$ID==mouse,]
  one.seq =one.seq[,!is.na(one.seq[1,])]
  maxtime = as.numeric(unlist(strsplit(colnames(one.seq)[ncol(one.seq)],"y"))[2])
  print(group_name)
  one.vec = generate_one(rt.days, parameter_vec=new_para_vec, maxtime=42, pd1=is.pd1group, group_name=old_names[group_name])
  #print(length(one.vec))
  one.round.simulated[one.round.simulated$ID==mouse,1:(maxtime)] = one.vec[1:maxtime]
  one.round.simulated$group[one.round.simulated$ID==mouse] = as.character(group_name)
  }
  

return(one.round.simulated)
}
