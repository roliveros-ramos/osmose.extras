
#' Create initialization file for an OSMOSE configuration
#'
#' @param input Filename of the main configuration file
#' @param output Ouput file containing the initialization configuration
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
initialize_osmose = function(input, output=NULL, ...) {
  
  conf = .readConfiguration(input)
  conf = .setupInitialization(conf)
  
  nsp = .getPar(conf, "simulation.nspecies")
  
  spind = .getPar(conf, "species.type") == "focal"
  spind = gsub(names(spind)[which(spind)], pattern="species.type.sp", replacement = "") 
  spnames = .getPar(conf, "species.name")[sprintf("species.name.sp%s", spind)]
  spind = as.numeric(spind)
  
  out = vector("list", nsp)
  names(out) = spnames
  pars = NULL
  
  for(sp in spind) {
    
    cat(sprintf("Initializing species %d\n", sp))
    this = .getPar(conf, sp=sp)
    sim = .simF_ini(conf, sp)
    sim$osmose = .initial_length_dist(sim, sp)
    pars = rbind(pars, as.matrix(sim$osmose))
    out[[.getPar(this, "species.name")]] = sim
  }
  
  pars = as.data.frame(pars)
  colnames(pars) = NULL
  
  xoutput = list(par=pars, init=out)
  class(xoutput) = c("osmose.initialization", class(xoutput))
  
  if(!is.null(output)) osmose::write_osmose(pars, file=output, sep=" = ")
    
  return(invisible(xoutput))
  
}

