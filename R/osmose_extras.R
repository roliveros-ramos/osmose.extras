
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
initialize_osmose = function(input, output=NULL, test=FALSE, ...) {
 
  ow = options("warn")
  options(warn=1)
  on.exit(options(ow))
  
  conf = .readConfiguration(input)
  # conf = .setupInitialization(conf)
  
  nsp = .getPar(conf, "simulation.nspecies")
  
  spind = .getPar(conf, "species.type") == "focal"
  spind = gsub(names(spind)[which(spind)], pattern="species.type.sp", replacement = "") 
  spnames = .getPar(conf, "species.name")[sprintf("species.name.sp%s", spind)]
  spind = as.numeric(spind)
  
  out = vector("list", nsp)
  names(out) = spnames
  pars = NULL
  
  for(sp in spind) {
    
    this = .getPar(conf, sp=sp)
    iSpName = .getPar(this, "species.name")
    
    cat(sprintf("\nInitializing species %d (%s)\n", sp, iSpName))
    
    sim = list()
    sim$cal       = read.cal(conf, sp)
    sim$biomass   = read.biomass(conf, sp)
    sim$yield     = read.yield(conf, sp)
    sim$fecundity = read.fecundity(conf, sp)
    sim$bioguess  = .getPar(.getPar(conf, sp=sp), "observed.biomass.guess")
    isp = sprintf("osmose.initialization.data.sp%d", sp)
    conf[[isp]]   = sim
    
    this = .getPar(conf, sp=sp)
    
    sim = .simF_ini(conf, sp, test=test)
    sim$osmose = .initial_length_dist(sim, sp)
    pars = rbind(pars, as.matrix(sim$osmose))
    # pars = rbind(pars, )
    out[[iSpName]] = sim
    
  }
  
  pars = as.data.frame(pars)
  colnames(pars) = NULL
  pars = pars[order(rownames(pars)), ]
  
  xoutput = list(par=pars, init=out)
  class(xoutput) = c("osmose.initialization", class(xoutput))
  
  if(!is.null(output)) osmose::write_osmose(pars, file=output, sep=" = ")
  
  return(invisible(xoutput))
  
}

