covinf.discrete = function(X, idx, response, dim = 1,
                           alpha, comp.dist, zero.prob, params.list)
{
  # Save a copy of X
  X.temp = X
  # Result to return
  df = data.frame(matrix(ncol=length(response),nrow=1+length(idx), 
                         dimnames=list(colnames(X[,c(1,idx)]), response)))
  
  # Base case
  X.temp[,idx] = 0
  for(j in 1:length(response)){
    df[1,j] = mean(
      switch (substr(response[j], 1, 2),
        "Me" = LRMoE::predict.mean(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
        "SD" = LRMoE::predict.var(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
        "VA" = LRMoE::predict.quantile(X.temp, alpha, comp.dist, zero.prob, params.list, 
                                       rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
        "CT" = LRMoE::predict.cte(X.temp, alpha, comp.dist, zero.prob, params.list, 
                                       rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
        # Error
        stop("Invalid input!")
      )
    )
  }
  
  # Other cases
  for(k in 1:length(idx))
  {
    X.temp[,idx[k]] = 1
    if(k>1){X.temp[,idx[k-1]] = 0}
    for(j in 1:length(response)){
      df[k+1,j] = mean(
        switch (substr(response[j], 1, 2),
                "Me" = LRMoE::predict.mean(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
                "SD" = LRMoE::predict.var(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
                "VA" = LRMoE::predict.quantile(X.temp, alpha, comp.dist, zero.prob, params.list, 
                                               rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
                "CT" = LRMoE::predict.cte(X.temp, alpha, comp.dist, zero.prob, params.list, 
                                          rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
                # Error
                stop("Invalid input!")
        )
      )
    }
    
  }
  
  return(df)
  
}



covinf.continuous = function(X, idx, eval.seq, response, dim = 1,
                           alpha, comp.dist, zero.prob, params.list)
{
  # Save a copy of X
  X.temp = X
  # Result to return
  df = data.frame(matrix(ncol=length(response),nrow=length(eval.seq), 
                         dimnames=list(paste(colnames(X)[idx], eval.seq, sep = ""), response)))
  
  # Run through eval.seq
  for(k in 1:length(eval.seq))
  {
    X.temp[,idx] = eval.seq[k]
    
    for(j in 1:length(response)){
      df[k,j] = mean(
        switch (substr(response[j], 1, 2),
                "Me" = LRMoE::predict.mean(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
                "SD" = LRMoE::predict.var(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
                "VA" = LRMoE::predict.quantile(X.temp, alpha, comp.dist, zero.prob, params.list, 
                                               rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
                "CT" = LRMoE::predict.cte(X.temp, alpha, comp.dist, zero.prob, params.list, 
                                          rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
                # Error
                stop("Invalid input!")
        )
      )
    }
    
  }
  
  return(df)
  
}













X.small = X[1:50,]
idx = c(15:20)
response = c("Mean", "VAR990", "CTE990")

# 
# df = data.frame(matrix(ncol=length(response),nrow=1+length(idx), 
#                        dimnames=list(colnames(X[,c(1,idx)]), response)))




temp = LRMoE::covinf.discrete(X.small, idx, response, dim = 1,
                       model.fit$alpha.fit, model.fit$comp.dist, model.fit$zero.fit, model.fit$params.fit)





X.small = X[1:50,]
idx = 2
eval.seq = seq(from = 0, to = 20, by = 1)
response = c("Mean", "VAR990", "CTE990")

temp = LRMoE::covinf.continuous(X.small, idx, eval.seq, response, dim = 1,
                         model.fit$alpha.fit, model.fit$comp.dist, model.fit$zero.fit, model.fit$params.fit)





