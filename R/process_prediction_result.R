#@ M: response probabilies, column names are the classes
#@ returns scores with class names
get_class_score = function(M){
  score = sapply(1:nrow(M), function(x){
          M[x,which.max(M[x,])];
      });
  score;
}

#@ M: response probabilies, column names are the classes
#@ family.list: list of methylation class families. Which family name as name, family member as vector 
get_family_score = function(M, family.list = MCF_family){
  #MCF = names(all_family)
  MCF_family_class = as.character(unlist(MCF_family))
  non_family = M[,which(!colnames(M) %in% MCF_family_class)]
  res = lapply(family.list, function(MCF){
    rowSums(M[,colnames(M) %in% MCF])
  })
  p = do.call(cbind, res)
  p = cbind(p, non_family)
  p
}

#@ x: vector
#@ family.list:
get_family_name = function(x, family.list){
  res = sapply(as.character(x), function(xx){
    if(xx %in% unlist(family.list)){
      mcf = lapply(family.list, function(MCF){ xx %in% MCF })
      names(mcf[which(mcf == 1)])
    }else{
      xx
    }
  })
  res
}
