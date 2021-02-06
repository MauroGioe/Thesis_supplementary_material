#Arguments
#Procedure: either "LASSO" or "Forward".
#Structure: either "Additive" or "Multiplicative".
#Inference: either "Controlled" for the conditional probability approach or "Uncontrolled" for carrying out inference without taking into account model selection.
#val.max: numeric value regarding the maximum value for an effect in the simulation. 
#The smallest effect is 0.1 and the other ones in the middle will assume increasing values adding up a quantity equal to (val.max - 0.1)/(n.var-1).
#n.var: numerical value that sets the number of true effects.
#n.fake: numerical value that sets the number of noise variables.
#n: numerical value that sets the number of observations.
#loops: numerical value that sets the number of loops for the simulation.
#noise.y: numerical value that sets the standard deviation of the random component of Y in the simulation.
#noise.x: numerical value that sets the standard deviation of the random component of X in the simulation.

#Functions to run before the simulations


check.function<-function(x)
{
  return(x<0.05)
}

tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),warning=W)
}


#Simulation
#There is no check for the validity of the arguments, be careful with the spelling

sim<-function(n=30,n.var=3,n.fake=3,loops=100,Procedure="LASSO",Structure="Additive",noise.y=1,noise.x=1,Inference="Uncontrolled",
                val.max=0.5)
{library(selectiveInference)
  library(glmnet)
  inter<-matrix(nrow=n,ncol=(factorial(n.var)/(factorial(n.var-2)*factorial(2) )  ) )
  pval.list<-list()
  pval.list.true<-list()
  pval.list.fake<-list()
  sel.alpha<-list()
  false.rej<-list()
  corr.rej<-list()
  Sig.increase<-list()
  TG.not<-0
  for(i in 1:loops )  {
    ### Identity matrix
    Sigma<-diag(n.var)
    Sigma.fake<-diag(n.fake)
    #Vector of means
    mu=rep(0,n.var)
    mu.fake=rep(0,n.fake)
    #Sampling
    x<-mvrnorm(n,mu=mu,Sigma=Sigma*noise.x)
    df <- data.frame(matrix(unlist(x), nrow=n, byrow=T),stringsAsFactors=FALSE)
    x.fake<-mvrnorm(n,mu=mu.fake,Sigma=Sigma.fake*noise.x)
    df.fake <- data.frame(matrix(unlist(x.fake), nrow=n, byrow=T),stringsAsFactors=FALSE)
    names(df.fake)<-paste(names(df.fake),".Fake")
    mu.predict<-rep(0,n.var+n.fake)
    Sigma.pred<-diag(n.var+n.fake)
    
    df.predict<-data.frame(matrix(unlist(mvrnorm(n,mu=mu.predict,Sigma=Sigma.pred*noise.x)), nrow=n, byrow=T),stringsAsFactors=FALSE)
    names(df.predict)[(n.var+1):(n.var+n.fake)]<-names(df.fake)
    df.fake<-cbind(df,df.fake)    
    #### Beta creation
    beta0<-5
    ###Y creation
    eps <- rnorm(n,0,sd=noise.y)
    nome<-rep("beta",n.var)
    nome<-paste(nome,seq(1,n.var,1))
    betas<-data.frame(t(seq(0.1,val.max,length.out =n.var)))
    betas.inte<-data.frame(t(seq(0.1,val.max,length.out = (factorial(n.var)/(factorial(n.var-2)*factorial(2) ) ) ) ) )
    names(betas)<-nome
    switch(Structure,"Additive"=
    {
      Y<-beta0+as.matrix(df)%*%(as.numeric(betas))+eps
      
    }, 
    "Multiplicative"={
      count<-0
      for (k in 1:(n.var-1)){
        for(j in (k+1):(n.var)){
          count<-count+1
          inter[,count]<-as.matrix(df)[,k]*as.matrix(df)[,j]
        }
      }
      Y<-beta0+as.matrix(df)%*%(as.numeric(betas))+inter%*%as.numeric(betas.inte)  +eps
    } )
    
    switch(Inference,"Controlled"= switch(Structure, 
                                         "Additive"= switch(Procedure,
                                                            "LASSO"={
                                                              mod<-model.matrix(~.,data=df.fake)
                                                              cv.fit <- cv.glmnet(mod[,-1], Y)
                                                              fit <- glmnet(mod[,-1], Y)
                                                              beta = coef(fit, s=cv.fit$lambda.1se/n, exact=TRUE,x=mod[,-1],y=Y)[-1]
                                                              out = fixedLassoInf(mod[,-1],Y,beta,cv.fit$lambda.1se, alpha = 0.05)
                                                              warn<-tryCatch.W.E(fixedLassoInf(mod[,-1],Y,beta,cv.fit$lambda.1se, alpha = 0.05))   
                                                              if(length(warn$warning)!=0){
                                                                if (grepl("TG.limits",warn$warning)==T) TG.not<-TG.not+1}
                                                              dati.pre<-model.matrix(~.,data=df.predict)
                                                              pval.list[[i]]<-data.frame(out$pv,names(out$vars))
                                                              pval.list[[i]]$ names.out.vars. <- gsub("`","",pval.list[[i]]$ names.out.vars.)
                                                              false.rej[[i]]<- sum(check.function(pval.list[[i]][grep(".Fake", pval.list[[i]][,2]),1]))
                                                              corr.rej[[i]]<- sum(!check.function(pval.list[[i]][grep(".Fake", pval.list[[i]][,2]),1]))                                                              
                                                              if (length(out$pv)<(n.var+n.fake)) {
                                                                pval.list[[i]]<-data.frame(out.pv=rep(1,length(names(df.fake))),names.out.vars.=names(df.fake))
                                                                pval.list[[i]]$out.pv[which( pval.list[[i]]$names.out.vars.%in% gsub("`","",names(out$vars)))]<-out$pv
                                                                
                                                              }
                                                            }
                                                            ,
                                                            "Forward"= {
                                                              mod<-model.matrix(~.,data=df.fake)
                                                              fsfit<-fs(mod[,-1],Y)
                                                              out.fs<-fsInf(fsfit,alpha = 0.05)
                                                              
                                                              warn<-tryCatch.W.E(fsInf(fsfit,alpha = 0.05))   
                                                              if(length(warn$warning)!=0){
                                                                if (grepl("TG.limits",warn$warning)==T) TG.not<-TG.not+1}
                                                              
                                                              
                                                              df.out<-as.data.frame(cbind(out.fs$vars,out.fs$pv))
                                                              colnames(df.out)<-c ("Vars", "pv")
                                                              ordered.out<-df.out[order(df.out$Vars),]
                                                              
                                                              dati.pre<-model.matrix(~.,data=df.predict)
                                                              pval.list[[i]]<-data.frame(out.pv=ordered.out$pv,names.out.vars.=colnames(mod)[-1])
                                                              pval.list[[i]]$ names.out.vars. <- gsub("`","",pval.list[[i]]$names.out.vars.)
                                                              false.rej[[i]]<- sum(check.function(pval.list[[i]][grep(".Fake", pval.list[[i]][,2]),1]))
                                                              corr.rej[[i]]<- sum(!check.function(pval.list[[i]][grep(".Fake", pval.list[[i]][,2]),1])) 
                                                              if (length(ordered.out$pv)<(n.var+n.fake)) {
                                                                pval.list[[i]]<-data.frame(out.pv=rep(1,length(names(df.fake))),names.out.vars.=names(df.fake))
                                                                pval.list[[i]]$out.pv[which( pval.list[[i]]$names.out.vars.%in% gsub("`","",names(df.fake[out.fs$vars])))]<-ordered.out$pv
                                                                
                                                              }
                                                            }
                                         )
                                         ,"Multiplicative"= switch(Procedure,
                                                                   "LASSO"={
                                                                     mod<-model.matrix(~.^2,data=df.fake)
                                                                     cv.fit <- cv.glmnet(mod[,-1], Y)
                                                                     fit <- glmnet(mod[,-1], Y)
                                                                     beta = coef(fit, s=cv.fit$lambda.1se/n, exact=TRUE,x=mod[,-1],y=Y)[-1]
                                                                     out = fixedLassoInf(mod[,-1],Y,beta,cv.fit$lambda.1se, alpha = 0.05)
                                                                     warn<-tryCatch.W.E(fixedLassoInf(mod[,-1],Y,beta,cv.fit$lambda.1se, alpha = 0.05))   
                                                                     if(length(warn$warning)!=0){
                                                                       if (grepl("TG.limits",warn$warning)==T) TG.not<-TG.not+1}
                                                                     
                                                                     
                                                                     df.out<-as.data.frame(cbind(out$vars,out$pv))
                                                                     colnames(df.out)<-c ("Vars", "pv")
                                                                     
                                                                     
                                                                     df.out.fake<-df.out[grep(".Fake",row.names(df.out)),]
                                                                     df.out.fake <-df.out.fake[-c(1:n.fake),]
                                                                     df.out.true<-df.out[-grep(".Fake",row.names(df.out)),]
                                                                     df.out.true<-df.out.true[-c(1:n.var),]
                                                                     
                                                                     dati.pre<-model.matrix(~.^2,data=df.predict)

                                                                     
                                                                     
                                                                     pval.list.true[[i]]<-data.frame(pv=df.out.true$pv)
                                                                     pval.list.fake[[i]]<-df.out.fake$pv
                                                                     false.rej[[i]]<- sum(check.function(pval.list.fake[[i]]))
                                                                     corr.rej[[i]]<- sum(!check.function(pval.list.fake[[i]]))
                                                                     
                                                                     if (length(df.out.true$pv)<length(betas.inte)){
                                                                       pval.list.true[[i]]<-data.frame(pv=rep(1,length(betas.inte)),vars=colnames(mod)[-grep(".Fake",colnames(mod))][-(1:(n.var+1))])
                                                                       pval.list.true[[i]]$pv[which( pval.list.true[[i]]$vars%in% gsub("`","",colnames(mod)[df.out.true$Vars+1]))]<-df.out.true$pv
                                                                     }
                                                                     if(length(df.out.fake$pv)< ( (factorial(n.var+n.fake)/(factorial(n.var+n.fake-2)*factorial(2) )  )-length(betas.inte) )    ){
                                                                       pval.list.fake[[i]]<-append(pval.list.fake[[i]],rep(1,(factorial(n.var+n.fake)/(factorial(n.var+n.fake-2)*factorial(2) )  )-length(betas.inte)-length(df.out.fake$pv) ))
                                                                       
                                                                     }
                                                                     
                                                                     
                                                                   }  
                                                                   ,
                                                                   "Forward"= {
                                                                     mod<-model.matrix(~.^2,data=df.fake)
                                                                     fsfit<-fs(mod[,-1],Y)
                                                                     out.fs<-fsInf(fsfit,alpha = 0.05)
                                                                     
                                                                     warn<-tryCatch.W.E(fsInf(fsfit,alpha = 0.05))   
                                                                     if(length(warn$warning)!=0){
                                                                       if (grepl("TG.limits",warn$warning)==T) TG.not<-TG.not+1}
                                                                     
                                                                     df.out<-as.data.frame(cbind(out.fs$vars,out.fs$pv))
                                                                     colnames(df.out)<-c ("Vars", "pv")
                                                                     df.out<-df.out[order(df.out$Vars),]
                                                                     df.out$Vars<-colnames(mod)[-1]
                                                                     
                                                                     dati.pre<-model.matrix(~.^2,data=df.predict)
                                                               
                                                                     
                                                                     df.out.fake<-df.out[grep(".Fake",df.out$Vars),]
                                                                     df.out.fake <-df.out.fake[-c(1:n.fake),]
                                                                     
                                                                     df.out.true<-df.out[-grep(".Fake",df.out$Vars),]
                                                                     df.out.true<-df.out.true[-c(1:n.var),]
                                                                     
                                                                     pval.list.true[[i]]<-data.frame(pv=df.out.true$pv,vars=df.out.true$Vars)
                                                                     pval.list.fake[[i]]<-df.out.fake$pv
                                                                     false.rej[[i]]<- sum(check.function(pval.list.fake[[i]]))
                                                                     corr.rej[[i]]<- sum(!check.function(pval.list.fake[[i]]))
                                                                     
                                                                     if (length(df.out.true$pv)<length(betas.inte)){
                                                                       pval.list.true[[i]]<-data.frame(pv=rep(1,length(betas.inte)),vars=colnames(mod)[-grep(".Fake",colnames(mod))][-(1:(n.var+1))])
                                                                       pval.list.true[[i]]$pv[which( pval.list.true[[i]]$vars%in% gsub("`","",colnames(mod)[df.out.true$Vars+1]))]<-df.out.true$pv
                                                                     }
                                                                     
                                                                     if(length(df.out.fake$pv)< ( (factorial(n.var+n.fake)/(factorial(n.var+n.fake-2)*factorial(2) )  )-length(betas.inte) )    ){
                                                                       pval.list.fake[[i]]<-append(pval.list.fake[[i]],rep(1,(factorial(n.var+n.fake)/(factorial(n.var+n.fake-2)*factorial(2) )  )-length(betas.inte)-length(df.out.fake$pv) ))
                                                                     }
                                                                     
                                                                     
                                                                   }
                                         )),
           "Uncontrolled"= switch(Structure,
                               "Additive"= switch(Procedure,
                                                  "LASSO"={
                                                    ols<-lm(Y~.,data=df.fake)
                                                    mod<-model.matrix(~.,data=df.fake)
                                                    cv.fit <- cv.glmnet(mod[,-1], Y)
                                                    fit <- glmnet(mod[,-1], Y)
                                                    Coefficients <- coef(fit, s = cv.fit$lambda.1se)
                                                    pval<-vector()
                                                    pval[which(tail(Coefficients,n=n.var+n.fake) == 0)]<-1
                                                    pval[which(tail(Coefficients,n=n.var+n.fake) != 0)]<-0
                                                    pval.list[[i]]<-data.frame(pval,rownames(Coefficients)[-1])
                                                    if (any(Coefficients[,1][ grep(".Fake",attributes(Coefficients[,1])$names)]!=0))
                                                      sel.alpha[[i]]<-1
                                                    else sel.alpha[[i]]<-0
                                                    Sig.increase[[i]]<- sum(Coefficients[,1][ grep(".Fake",attributes(Coefficients[,1])$names)]!=0)-sum((summary(ols)$coefficients[,4])[grep(".Fake",attributes(summary(ols)$coefficients[,4])$names)]<0.05)
                                                  }
                                                  ,
                                                  "Forward"={
                                                    mod<-(lm(Y~.,data=df.fake ))
                                                    min.mod<-lm(Y~1,data=df.fake)
                                                    st<-summary(step(min.mod,direction = "forward",scope=formula(mod),trace = 0))
                                                    sel.alpha[[i]]<- sum( st$coefficients[,4][ grep(".Fake",attributes(st$coefficients[,1])$names)]<0.05)/length(st$coefficients[,4][ grep(".Fake",attributes(st$coefficients[,1])$names)])
                                                    if(sel.alpha[[i]]=="NaN") sel.alpha[[i]]<-0
                                                    df.pval<-data.frame(Coef=seq(1:ncol(df.fake)),Pval=rep(1,ncol(df.fake)))
                                                    Sig.increase[[i]]<- sum( st$coefficients[,4][ grep(".Fake",attributes(st$coefficients[,1])$names)]<0.05)- sum((summary(mod)$coefficients[,4])[grep(".Fake",attributes(summary(mod)$coefficients[,4])$names)]<0.05)
                                                    if(st$df[1]>1){
                                                      obj<- row.names(st$coefficients)[2:length(row.names(st$coefficients))]
                                                      obj.true<-obj[-grep(".Fake",obj)]
                                                      obj.fake<-obj[grep(".Fake",obj)]
                                                      
                                                      matches.true <- as.numeric(regmatches(obj.true, gregexpr("\\d+", obj.true)))
                                                      matches.fake <- as.numeric(regmatches(obj.fake, gregexpr("\\d+", obj.fake)))
                                                      
                                                      matches<-c(matches.true,(matches.fake+n.var))
                                                      
                                                      min.p<-2+(st$df[1]*3)
                                                      max.p<-st$df[1]+(st$df[1]*3)
                                                      
                                                      df.pval$Pval[matches]<- st$coefficients[c(min.p:max.p)]
                                                      pval.list[[i]]<-data.frame(pval=df.pval$Pval,vars=seq(1:ncol(df.fake)))
                                                      
                                                    }
                                                    else pval.list[[i]]<-data.frame(pval=df.pval$Pval,vars=seq(1:ncol(df.fake)))
                                                  }
                               )
                               ,"Multiplicative"= switch(Procedure,
                                                         "LASSO"={
                                                           ols<-lm(Y~.^2,data=df.fake)
                                                           mod<-model.matrix(~.^2,data=df.fake)
                                                           cv.fit <- cv.glmnet(mod[,-1], Y)
                                                           fit <- glmnet(mod[,-1], Y)
                                                           Coefficients <- coef(fit, s = cv.fit$lambda.1se)
                                                           
                                                           df.out.fake<-Coefficients[grep(".Fake",row.names(Coefficients)),]
                                                           df.out.fake <-df.out.fake[-c(1:n.fake)]
                                                           df.out.fake<-ifelse(df.out.fake==0,1,0)
                                                           df.out.true<-Coefficients[-grep(".Fake",row.names(Coefficients)),]
                                                           df.out.true<-df.out.true[-c(1:(n.var+1))]
                                                           df.out.true<-ifelse(df.out.true==0,1,0)
                                                           
                                                           if (any(df.out.fake==0))
                                                             sel.alpha[[i]]<-1
                                                           else sel.alpha[[i]]<-0
                                                           Sig.increase[[i]]<- sum(df.out.fake==0)-sum((summary(ols)$coefficients[,4])[grep(".Fake",attributes(summary(ols)$coefficients[,4])$names)]<0.05)
                                                           
                                                           pval.list.true[[i]]<-data.frame(pv=df.out.true)
                                                           pval.list.fake[[i]]<-df.out.fake
                                                           
                                                           
                                                         }
                                                         
                                                         ,
                                                         "Forward"={
                                                           mod<-(lm(Y~.^2,data=df.fake ))
                                                           min.mod<-lm(Y~.,data=df.fake)
                                                           st<-summary(step(min.mod,direction = "forward",scope=formula(mod),trace = 0))
                                                           sel.alpha[[i]]<-sum(check.function(st$coefficients[,4][grepl("(?=.*.Fake)(?=.*:)",attributes(st$coefficients[,1])$names,perl = T)]))/length(st$coefficients[,4][grepl("(?=.*.Fake)(?=.*:)",attributes(st$coefficients[,1])$names,perl = T)])
                                                           if(sel.alpha[[i]]=="NaN") sel.alpha[[i]]<-0
                                                           
                                                           n.inter<-length(mod$coefficients)-1-n.var-n.fake
                                                           df.pval<-data.frame(Coef=attributes(mod$coefficients)$names[(2+n.var+n.fake):(n.inter+1+n.var+n.fake)],Pval=rep(1,n.inter))
                                                           
                                                           
                                                           
                                                           df.out.fake<-df.pval[grep(".Fake",df.pval$Coef),]
                                                           df.out.fake$Pval[which(df.out.fake$Coef%in%attributes(st$coefficients[,1])$names )]<-st$coefficients[,4][which(attributes(st$coefficients[,1])$names %in%   df.out.fake$Coef)]
                                                           df.out.true<-df.pval[-grep(".Fake",df.pval$Coef),]
                                                           df.out.true$Pval[which(df.out.true$Coef%in%attributes(st$coefficients[,1])$names )]<-sort(st$coefficients[,4][which(attributes(st$coefficients[,1])$names %in%   df.out.true$Coef)],T)
                                                           
                                                           pval.list.true[[i]]<-data.frame(pv=df.out.true$Pval)
                                                           pval.list.fake[[i]]<-df.out.fake$Pval
                                                           Sig.increase[[i]]<- sum(check.function(st$coefficients[,4][grepl("(?=.*.Fake)(?=.*:)",attributes(st$coefficients[,1])$names,perl = T)])) -sum((summary(mod)$coefficients[,4])[grep(".Fake",attributes(summary(mod)$coefficients[,4])$names)]<0.05)
                                                           
                                                         }
                               ))
    )
    
  }      
  
  switch(Structure,"Additive"={ 
    pval.list<-lapply(pval.list, function(x) x[-2])
    results<-matrix(unlist(pval.list),nrow=n.var+n.fake)},
    "Multiplicative"={
      
      results.pow<-matrix(unlist(pval.list.true),nrow=length(betas.inte))
      
      
      if(Inference=="Uncontrolled"& Procedure=="Forward")
        results.alpha<-matrix(unlist(pval.list.fake),nrow= length(mod$coefficients)-1-n.var-n.fake-length(betas.inte)    )
      else
        results.alpha<-matrix(unlist(pval.list.fake),nrow= ncol(mod[,-1])-n.var-n.fake-length(betas.inte)    )
    } )
  
  switch(Structure,"Additive"={
    pow<-apply(apply(results[1:n.var,],2,check.function),1,sum)/loops
    if (Inference=="Uncontrolled") {alpha<-apply(apply(results[(n.var+1):(n.var+n.fake),],2,check.function),1,sum)/loops
    alpha<-mean(alpha)
    }
  }
  ,"Multiplicative"={
    
    pow<-apply(apply(results.pow,2,check.function),1,sum)/loops
  }
  )
  
  if (Inference=="Uncontrolled") Sig.increase<-mean(unlist(Sig.increase))
  pow<-as.data.frame(t(pow))
  if (Inference=="Controlled")  sel.alpha<-sum(unlist(false.rej))/((sum(unlist(corr.rej)))+sum(unlist(false.rej)))
  
  switch(Structure,"Additive"={colnames(pow)<-betas},"Multiplicative"={colnames(pow)<-betas.inte})
  if (Inference=="Controlled"){ Res<-list(Power=pow,sel.alpha=mean(unlist(sel.alpha)),TG.not= (TG.not/loops)*100  )}
  else{
  Res<-list(Power=pow,sel.alpha=mean(unlist(sel.alpha)),Sig.increase=Sig.increase)}
  return(Res)
  
}

