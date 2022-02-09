
surveyBiomass = function(object, conf) {
  
  x = get_var(object, "biomassBySize")
  
  size = as.numeric(colnames(x[[1]]))
  thr = .getPar(conf, par="biomass.cutoff.size")
  nsp = .getPar(conf, "simulation.nspecies")
  spind = .getPar(conf, "species.type") == "focal"
  spind = gsub(names(spind)[which(spind)], pattern="species.type.sp", replacement = "") 
  spnames = .getPar(conf, "species.name")[sprintf("species.name.sp%s", spind)]
  spind = setNames(as.numeric(spind), nm=spnames)
  
  nrow = nrow(x[[1]])
  ncol = length(x)
  nrep = ifelse(is.na(dim(x[[1]])[3]), 1, dim(x[[1]][3]))
  
  out = array(dim=c(nrow, ncol, nrep))
  
  for(i in seq_along(x)) {
    
    isp = spind[names(x)[i]]
    ix = x[[i]]
    if(nrep==1) dim(ix) = c(dim(ix), 1)
    ind = size >= .getPar(thr, sp=isp)
    out[, i, ] = apply(ix[, ind, , drop=FALSE], c(1,3), FUN=sum)
    
  }
  
  rownames(out) = rownames(x[[1]])
  colnames(out) = names(x)
  
  class(out) = c("osmose.biomass", "array")
  return(out)
  
}
