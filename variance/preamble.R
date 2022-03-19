

library(data.table)


##code to generate one mouse given a set of parameters
#source('~/Documents/Dissertation/PAI/app/dyn_sys_equations.R')
#source('~/Documents/Dissertation/PAI/app/verification_code.R')

# totaltime<-40
# setwd("/home/marylena/Documents/Dissertation/PAI/variance_determination")
# source("fit_one_dataset.R")
# source("generate_one_dataset.R")
# source("generate_one.R")
one.round = read.csv("/home/marylena/Documents/Dissertation/data/log-raw-data/Round1_log.csv")[,-c(2)]
colnames(one.round)[1]="ID"

one.group = one.round[one.round$group=="0Gy_PDL"]


days = strsplit(colnames(one.round),"y")
days = as.numeric(unlist(days))[!is.na(as.numeric(unlist(days)))]

#noise_vec = c(0.002, 0.16,0.1,0.001,0.002,0.001,0.001, 0.2)
#IC_noise = rep(NA, length(names(Tn_loop)))
#names(IC_noise) = names(Tn_loop)
#names(noise_vec) = c("mu","rho","lambda","omega1","omega2","BT0","TT0", "Tinit")

old_names= levels(one.round$group)
# new_names = names(Tn_loop)[c(1:8,11,12,9,10)]
# names(old_names) = new_names
# names(new_names)=old_names
# real_params = fit_one_dataset(one.round, splitchar = "y")

# starting_params = noise_vec[2:7]
# penalty_tracker=c()
# obj_fnc = function(noise_vec, mu){
#   noise_vec_01 = c(mu,abs(noise_vec))
#   print(noise_vec_01)
#   simulated_dataset=generate_one_dataset(noise_vec_01)
#   simulated_params = fit_one_dataset(simulated_dataset)
#   penalty = sum((simulated_params[1:2]-real_params[1:2])^2+900000*(simulated_params[2]=0))
#   penalty_tracker=c(penalty_tracker, penalty)
#   return(penalty)
# }
# 
# optim_result = optim(par=starting_params,fn=obj_fnc, mu=.002, method="L-BFGS-B")
# simulated_dataset=generate_one_dataset(c(.002,optim_result$par))
# simulated_params = fit_one_dataset(simulated_dataset)
