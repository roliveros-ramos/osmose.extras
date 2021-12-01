
# Simulation --------------------------------------------------------------


.simCatch_ini = function(conf, sp) {
  
  this = .getPar(conf, sp=sp)
  ndt = conf$simulation.time.ndtperyear
  A = .getPar(this, "species.lifespan")
  a = .getPar(this, "species.length2weight.condition.factor")
  b = .getPar(this, "species.length2weight.allometric.power")
  
  age_bins = seq(from=0, to=A, by=1/ndt)
  age = 0.5*head(age_bins, -1) + 0.5*tail(age_bins, -1)
  size = VB(age, this, method=3)
  size_bins = VB(age_bins, this, method=3)
  
  T = 5*length(age)
  C = length(age)
  
  CAL = read.cal(conf, sp=sp)
  if(is.null(CAL)) return(NULL)
  harvested = CAL$harvested
  BIO = read.biomass(conf, sp)
  fecundity = read.fecundity(conf, sp)
  
  rF = fecundity[1:ndt]/sum(fecundity[1:ndt])
  bioguess = .bioguess(x=BIO, ndt=ndt)
  
  
  
  trans = rebinning(CAL$bins, VB(age_bins, this, method=3))
  cal = CAL$mat %*% trans
  weight = a*size^b
  yield_obs   = rep((1e-6*rowSums(t(t(cal)*weight)))[1:ndt], length=T)
  
  M = calculateMortality(conf, sp)
  Ma = M$M[cut(age, breaks = M$age, labels = FALSE)]
  
  
  xn = which.max(size >= .getPar(this, "observed.biomass.cutoff.size")) - 1
  
  ini = exp(-cumsum(c(0,Ma[-length(Ma)]/ndt)))
  inibio = 1e-6*sum(ini*weight)
  R = ndt*bioguess/inibio
  ini = R*ini/sum(ini)
  
  if(bioguess==0) {
    # species not starting in the simulation
    output = list(pop=rep(0, C), catch=rep(0, C), R=0, 
                  biomass=rep(0, ndt), yield=rep(0, ndt), F=rep(0, ndt),
                  dist = rep(0, C), distB = rep(0, C), 
                  selectivity=rep(0, C), age=age, size=size,
                  Fguess=0, observed=list(biomass=bioguess, yield=yield_obs[1:ndt]),
                  bins=list(age=age_bins, size=size_bins), harvested=FALSE)
    
    class(output) = c("osmose.init",class(output))
    
    return(output)
    
  }
  
  for(ix in 1:5) {
    
    pop = matrix(0, nrow=T, ncol=C)
    catch = matrix(NA, nrow=T, ncol=C)
    pop[1, ] = ini
    pop[ , 1] = R*rep(rF, length=T)
    
    for(t in 1:T) {
      it = ((t-1) %% ndt) + 1
      MID = pop[t, ]*exp(-0.5*Ma/ndt) # first half
      MID.C = pmax(MID - cal[it, ], 0.01*MID) # remove as much catch as possible, CHECK
      catch[t, ] = MID - MID.C
      END = MID.C*exp(-0.5*Ma/ndt) # second half
      if(t!=T) pop[t+1, 2:C] = head(END, -1)   
    }
    
    ixn = seq_len(ncol(pop))
    ixn = if(xn>0) tail(ixn, -xn) else ixn 
    biopop = 1e-6*t(t(pop)*weight)[, ixn]
    biomass = rowSums(biopop)
    bioref = mean(tail(biomass, ndt))
    yield   = 1e-6*rowSums(t(t(catch)*weight))
    
    Fguess = tail(yield_obs, ndt)/tail(biomass, ndt)
    Fseason = Fguess/sum(Fguess)
    if(all(is.na(Fseason))) Fseason = rep(1, ndt)/ndt
    
    Rupb = bioguess/bioref
    Rupc = mean(tail(yield_obs/yield, ndt))
    
    inc = weighted.mean(c(Rupb, Rupc), w = c(1, 2), na.rm=TRUE)
    R = R*inc
    
    ini = inc*pop[T - ndt + 1 ,]
    
  }
  
  xind = T - ndt + 1
  pind = xind + seq_len(ndt) - 1
  
  es = colSums(catch[pind, ])/colMeans(pop[pind, ])
  
  Fguess = sum(yield[pind]/biomass[pind])
  
  output = list(pop=pop[pind, ], catch=catch[pind, ], R=R, 
                biomass=biomass[pind], yield=yield[pind], F=Fseason,
                dist = pop[xind ,], distB = 1e-6*pop[xind ,]*weight, 
                selectivity=es, age=age, size=size,
                Fguess=Fguess, observed=list(biomass=bioguess, yield=yield_obs[1:ndt]),
                bins=list(age=age_bins, size=size_bins), harvested=harvested)
  
  class(output) = c("osmose.init",class(output))
  return(output)
  
}
.simF_ini = function(conf, sp, tiny=1e-3, cv=c(0.1, 0.1)) {
  
  this = .getPar(conf, sp=sp)
  ndt = conf$simulation.time.ndtperyear
  A = .getPar(this, "species.lifespan")
  a = .getPar(this, "species.length2weight.condition.factor")
  b = .getPar(this, "species.length2weight.allometric.power")
  
  age_bins = seq(from=0, to=A, by=1/ndt)
  age = 0.5*head(age_bins, -1) + 0.5*tail(age_bins, -1)
  size = VB(age, this, method=3)
  size_bins = VB(age_bins, this, method=3)
  
  T = 5*length(age)
  C = length(age)
  
  CAL = read.cal(conf, sp=sp)
  harvested = if(is.null(CAL)) TRUE else CAL$harvested
  if(!isTRUE(harvested)) return(.simCatch_ini(conf, sp))
  BIO = read.biomass(conf, sp)
  fecundity = read.fecundity(conf, sp)
  
  rF = fecundity[1:ndt]/sum(fecundity[1:ndt])
  biofit   = .bioguess(x=BIO, ndt=ndt, ts=TRUE)
  bioguess = .bioguess(x=BIO, ndt=ndt)
  
  yield_obs   = read.yield(conf, sp)
  yield_obs   = rep(yield_obs[1:ndt], length=T)
  
  weight = a*size^b
  M = calculateMortality(conf, sp)
  Ma = M$M[cut(age, breaks = M$age, labels = FALSE)]
  
  xn = which.max(size >= .getPar(this, "observed.biomass.cutoff.size")) - 1
  
  sim = .simCatch_ini(conf, sp)
  
  if(is.null(sim)) {
    
    sel = .getSelectivity(size, this) # creado con parametros
    Fseason = yield_obs[1:ndt]/sum(yield_obs[1:ndt])
    
    ini = exp(-cumsum(c(0,Ma[-length(Ma)]/ndt)))
    inibio = 1e-6*sum(ini*weight)
    R = ndt*bioguess/inibio
    xdist = R*ini/sum(ini)
    
    sim = list(R=R, Fguess=sum(yield_obs[1:ndt]/bioguess))
    
    
  } else {
    
    es = empirical_selectivity(matrix(sim$selectivity, nrow=1), fleet = "sim",
                               years = 1, bins = sim$size)
    ss_emp = suppressMessages(fit_selectivity(es, pattern=27, k=5))
    sel = pmax(c(ss_emp$selectivity), tiny)
    Fseason = sim$F
    xdist = sim$dist
  }
  
  .simF = function(par, value=FALSE) {
    
    R = exp(par[1])
    F = exp(par[2])
    
    pop = matrix(0, nrow=T, ncol=C)
    catch = matrix(NA, nrow=T, ncol=C)
    pop[1, ] = xdist
    pop[ , 1] = R*rep(rF, length=T)
    
    for(t in 1:T) {
      it = ((t-1) %% ndt) + 1
      Ft = F*Fseason[it]*sel
      Zt = Ma/ndt + Ft
      tmp = pop[t, ]*exp(-Zt) 
      catch[t, ] = pop[t, ]*(Ft/Zt)*(1-exp(-Zt)) 
      if(t!=T) pop[t+1, 2:C] = head(tmp, -1)   
    }
    
    ixn = seq_len(ncol(pop))
    ixn = if(xn>0) tail(ixn, -xn) else ixn 
    biopop = 1e-6*t(t(pop)*weight)[, ixn]
    biomass = rowSums(biopop)
    yield   = 1e-6*rowSums(t(t(catch)*weight))
    
    ll_biomass  = lnorm2(biofit, tail(biomass, ndt))
    ll_yield    = lnorm2(yield_obs[1:ndt] , tail(yield, ndt))
    
    ll = ll_biomass*llw(cv[1]) + ll_yield*llw(cv[2])
    
    if(!isTRUE(value)) return(ll)
    
    xind = T - ndt + 1
    pind = xind + seq_len(ndt) - 1
    
    es = colSums(catch[pind, ])/colMeans(pop[pind, ])
    
    output = list(pop=pop[pind, ], catch=catch[pind, ], R=R, 
                  biomass=biomass[pind], yield=yield[pind], F=Fseason, Fguess=F,
                  dist=pop[xind ,], distB = 1e-6*pop[xind, ]*weight, 
                  selectivity=sel, age=age, size=size,
                  observed=list(biomass=biofit, yield=yield_obs[1:ndt]),
                  bins=list(age=age_bins, size=size_bins))
    
    return(output)
    
  }
  
  opt = calibrar::calibrate(par=log(c(sim$R, sim$Fguess)), 
                            fn = .simF, method = "L-BFGS-B")
  
  output = c(.simF(opt$par, value=TRUE),  opt=list(opt))
  class(output) = c("osmose.init",class(output))
  return(output)
  
}


# Internal ----------------------------------------------------------------


lnorm2 = function(obs, sim, tiny=1e-2, ...) {
  if(all(!is.finite(sim))) return(Inf)
  obs = log(obs + tiny)
  sim = log(sim + tiny)
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

llw = function(cv) 1/(2*cv^2)

.initial_length_dist = function(sim, sp) {
  
  dist = sim$distB
  bio_ini = sum(dist)
  bio_rel = if(bio_ini==0) dist else dist/bio_ini
  bins = sim$bins$size
  tl_sp = rep(2, length(dist))
  
  out = c(round(bio_ini, 1), 
          paste(format(bio_rel, scientific = FALSE), collapse=", "), 
          paste(round(head(bins, -1), 2), collapse=","),
          paste(round(tl_sp, 2), collapse=", "))
  dim(out) = c(4, 1)
  
  out = as.data.frame(out)
  rownames(out) = sprintf(c("population.initialization.biomass.sp%d",
                            "population.initialization.relativebiomass.sp%d",
                            "population.initialization.size.sp%d",
                            "population.initialization.tl.sp%d"), sp)
  colnames(out) = NULL
  
  return(out)
  
}

