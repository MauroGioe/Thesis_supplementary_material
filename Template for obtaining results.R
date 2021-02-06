library(abind)
#Change the parameters of the sim function according to your needs
#The warnings for Inference="Controlled" are related to the TG constraints not being satisfied

# Template - Results for Inference="Uncontrolled"
n.loops<-c(1000,1000,1000,1000)
list_results<-list()
results_pow<-list()
results_sel_alpha<-list()
results_sig_increase<-list()

for(i in 1:length(n.loops)){
  list_results[[i]]<-sim(loops = n.loops[i],Inference = "Uncontrolled")
    results_pow[[i]]<-list_results[[i]]$Power
    results_sel_alpha[[i]]<-list_results[[i]]$sel.alpha
    results_sig_increase[[i]]<-list_results[[i]]$Sig.increase
print(paste("Loop",i,"completed"))
    }

#mean
result_pow_mean<-apply(do.call(mapply, c('abind', results_pow, rev.along = 0)),2,mean)
result_sel_alpha_mean<-apply(do.call(mapply, c('abind', results_sel_alpha, rev.along = 0)),2,mean)
result_sig_increase_mean<-apply(do.call(mapply, c('abind', results_sig_increase, rev.along = 0)),2,mean)
#range
result_pow_range<-apply(do.call(mapply, c('abind', results_pow, rev.along = 0)),2,range)
result_sel_alpha_range<-apply(do.call(mapply, c('abind', results_sel_alpha, rev.along = 0)),2,range)
result_sig_increase_range<-apply(do.call(mapply, c('abind', results_sig_increase, rev.along = 0)),2,range)


# Template - Results for Inference="Controlled"

n.loops<-c(1000,1000,1000,1000)
list_results<-list()
results_pow<-list()
results_sel_alpha<-list()
results_TG_constr<-list()

for(i in 1:length(n.loops)){
  list_results[[i]]<-sim(loops = n.loops[i],Inference = "Controlled")
  results_pow[[i]]<-list_results[[i]]$Power
  results_sel_alpha[[i]]<-list_results[[i]]$sel.alpha
  results_TG_constr[[i]]<-list_results[[i]]$TG.not
  print(paste("Loop",i,"completed"))
}

#mean
result_pow_mean<-apply(do.call(mapply, c('abind', results_pow, rev.along = 0)),2,mean)
result_sel_alpha_mean<-apply(do.call(mapply, c('abind', results_sel_alpha, rev.along = 0)),2,mean)
results_TG_constr_mean<-apply(do.call(mapply, c('abind', results_TG_constr, rev.along = 0)),2,mean)

#range
result_pow_range<-apply(do.call(mapply, c('abind', results_pow, rev.along = 0)),2,range)
result_sel_alpha_range<-apply(do.call(mapply, c('abind', results_sel_alpha, rev.along = 0)),2,range)
results_TG_const_range<-apply(do.call(mapply, c('abind', results_TG_constr, rev.along = 0)),2,range)
