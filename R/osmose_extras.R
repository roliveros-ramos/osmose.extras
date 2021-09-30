
.readConfiguration = function(file, usePath=TRUE, ...) {
  
  .guessSeparator = function(Line){
    SEPARATORS = c(equal = "=", semicolon = ";",
                   coma = ",", colon = ":", tab = "\t")
    guess = which.min(nchar(lapply(str_split(Line,SEPARATORS), "[", i = 1)))
    separator = SEPARATORS[guess]
    
    return(separator)
  }
  
  .getKey = function(Line, KeySeparator) {
    Key = str_split(Line, KeySeparator)[[1]][1]
    return(stringr::str_trim(Key))
  }
  
  .getValues = function(x, KeySeparator){
    start = str_locate(x, pattern=KeySeparator)[1,1]
    if(is.na(start)) return(NULL)
    values = stringr::str_sub(x, start+1, nchar(x))
    valueseparator = .guessSeparator(values)
    values = stringr::str_trim(str_split(values, valueseparator)[[1]])
    values = values[nchar(values)!=0]
    values = .guessType(values)
    return(values)
  }
  
  .guessType = function(x) {
    xx = tolower(x)
    if(identical(xx, "null")) return(NULL)
    if(identical(xx, "NA")) return(NULL)
    x[xx=="true"] = "TRUE"
    x[xx=="false"] = "FALSE"
    x = x[xx!="null"]
    out = type.convert(c(x), as.is = TRUE)
    return(out)
  }
  
  .comment_trim = function(x, char="#") {
    start = str_locate(x, pattern=char)[1,1]
    if(is.na(start)) return(x)
    return(str_sub(x, 1, start - 1))
  }
  
  .addPath = function(x, path) {
    if(is.null(x)) return(x)
    if(!is.character(x)) return(x)
    # if(!is.null(attr(x, "path"))) 
      # path = file.path(attr(x, "path"), path)
    if(file.exists(file.path(path, x))) 
      attr(x, "path") = normalizePath(path, winslash = "/", mustWork = FALSE)
    return(x)
  }
  
  if(!is.null(attr(file, "path"))) file = c(file.path(attr(file, "path"), file))

  config = readLines(file) # read lines
  config = lapply(config, .comment_trim) # remove comments
  config = lapply(config, str_trim)
  config[grep("^[[:punct:]]", config)] = NULL
  config = config[nchar(config)!=0]
  
  keySeparator  = sapply(config, .guessSeparator)
  key           = mapply(.getKey, config, keySeparator)
  values        = mapply(.getValues, config, keySeparator, SIMPLIFY = FALSE)
  
  names(values) = tolower(key)
  
  if(isTRUE(usePath)) values = lapply(values, .addPath, path = dirname(file))
  
  ii = grep(names(values), pattern="osmose.configuration")
  if(length(ii) == 0) return(values)
  
  while(length(ii) > 0) {
    
    values = append(values, .readConfiguration(values[[ii[1]]]), ii[1])
    values = values[-ii[1]]
    ii = grep(names(values), pattern="osmose.configuration")
  }
 
  return(values)
}

.getPar = function(conf, par, sp=NULL) {
  if(!is.null(sp)) par = sprintf(".sp%d$", sp)
  par = tolower(par)
  out = conf[grep(names(conf), pattern=par)]
  if(length(out)==1) out = out[[1]]
  if(length(out)==0) out = NULL
  return(out)
}


VB = function(age, this, method=3) {
  
  k = .getPar(this, "species.k") 
  Linf = .getPar(this, "species.linf")
  t0 = .getPar(this, "species.t0")
  t_thr = .getPar(this, "species.vonbertalanffy.threshold.age")
  L_egg = .getPar(this, "species.egg.size")
  L_thr = Linf*(1 - exp(-k*(t_thr - t0)))
  
  if(method==0) {
    # original VB
    l = ifelse(age < 0, NA, Linf*(1 - exp(-k*(age - t0))))
    return(l)
  }
  if(method==1) {
    # linear approximation (current OSMOSE)
    a = (L_thr - L_egg)/t_thr
    b = L_egg
    l = ifelse(age < t_thr, a*age + b, Linf*(1 - exp(-k*(age - t0))))
    return(l)
  }
  if(method==2) {
    # cuadratic (match up to first derivative)
    a = k*(Linf - L_thr)/t_thr - L_thr/(t_thr)^2
    b = -k*(Linf - L_thr) + 2*L_thr/t_thr
    c = L_egg
    l = ifelse(age < t_thr, a*age^2 + b*age + c, Linf*(1 - exp(-k*(age - t0))))
    return(l)
  }
  if(method==3) {
    # cubic (match up to second derivative)
    a = (L_thr - L_egg)/t_thr^3 - k*(Linf - L_thr)/t_thr^2 - (k^2/2)*(Linf - L_thr)/t_thr
    b = (k^2 + 3*k/t_thr)*(Linf - L_thr) - 3*(L_thr - L_egg)/t_thr^2  
    c = -2*k*(Linf - L_thr) - (k^2/2)*(Linf - L_thr)*t_thr + 3*(L_thr - L_egg)/t_thr
    d = L_egg
    l = ifelse(age < t_thr, a*age^3 + b*age^2 + c*age + d, Linf*(1 - exp(-k*(age - t0))))
    return(l)
  }
  stop("Incorrect 'method' specification.")
}

VB_inv = function(size, this, method=3) {
  
  k = .getPar(this, "species.k") 
  Linf = .getPar(this, "species.linf")
  t0 = .getPar(this, "species.t0")
  t_thr = .getPar(this, "species.vonbertalanffy.threshold.age")
  L_egg = .getPar(this, "species.egg.size")
  L_thr = Linf*(1 - exp(-k*(t_thr - t0)))
  A = .getPar(this, "species.lifespan")
  
 .invVB = function(size, k, Linf, t0, A) {
    out = suppressWarnings((-1/k)*log(1-size/Linf) + t0)
    out[out>A | is.na(out)] = A
    return(out)
  }
  
  if(method==0) {
    # original VB
    age = ifelse(size < 0, NA, .invVB(size, k, Linf, t0, A))
    return(age)
  }
  if(method==1) {
    # linear approximation (current OSMOSE)
    a = (L_thr - L_egg)/t_thr
    b = L_egg
    age = ifelse(size < L_thr, (size-b)/a, .invVB(size, k, Linf, t0, A))
    return(age)
  }
  if(method==2) {
    # cuadratic (match up to first derivative)
    a = k*(Linf - L_thr)/t_thr - L_thr/(t_thr)^2
    b = -k*(Linf - L_thr) + 2*L_thr/t_thr
    c = L_egg - size
    athr = suppressWarnings(cbind((-b + sqrt(b^2 - 4*a*c))/(2*a),
                 (-b - sqrt(b^2 - 4*a*c))/(2*a)))
    athr[athr>t_thr | is.nan(athr)] = NA
    athr = suppressWarnings(apply(athr, 1, min, na.rm=TRUE))
    age = ifelse(size < L_thr, athr, .invVB(size, k, Linf, t0, A))
    return(age)
  }
  if(method==3) {
    # cubic (match up to second derivative)
    a = (L_thr - L_egg)/t_thr^3 - k*(Linf - L_thr)/t_thr^2 - (k^2/2)*(Linf - L_thr)/t_thr
    b = (k^2 + 3*k/t_thr)*(Linf - L_thr) - 3*(L_thr - L_egg)/t_thr^2  
    c = -2*k*(Linf - L_thr) - (k^2/2)*(Linf - L_thr)*t_thr + 3*(L_thr - L_egg)/t_thr
    d = L_egg - size
    # solution to the cubic, guessing is the real one!
    D0 = b^2 - 3*a*c
    D1 = 2*b^3 - 9*a*b*c + 27*a^2*d
    C = ((D1 + sqrt(D1^2 - 4*D0^3))/2)^(1/3)
    athr = -(b + C + D0/C)/(3*a)
    age = ifelse(size < L_thr, athr, .invVB(size, k, Linf, t0, A))
    return(age)
  }
  stop("Incorrect 'method' specification.")
}

calculateMLF = function(conf, sp) {
  
  .calculateMLF = function(fecundity, isMature, weight) {
    
    .mlf = function(delay, fecundity, isMature, weight) {
      MLF = sum(isMature*fecundity[delay + seq_len(nlife)]*weight)
      return(MLF)
    }
    
    nfec  = length(fecundity)
    nlife = length(weight)
    if(nfec <= nlife) {
      fecundity = rep(fecundity, length=nlife+nfec)
      niter = nfec
    } else {
      niter = nfec - nlife
    }
    nfec = length(fecundity)
    
    out = sapply(0:niter, .mlf, fecundity=fecundity, isMature=isMature, weight=weight)
    return(out)
  }
  
  ndt = conf$simulation.time.ndtperyear
  this = .getPar(conf, sp=sp)
  
  a = .getPar(this, "species.length2weight.condition.factor")
  b = .getPar(this, "species.length2weight.allometric.power")
  tn = .getPar(this, "species.lifespan")
  
  age = seq(from=0+0.5/ndt, to=tn, by=1/ndt)
  size = VB(age, this, method=3)
  weight = a*size^b
  
  repfile = .getPar(this, "reproduction.season.file")
  fecundity = as.numeric(unlist(read.csv(file.path(attr(repfile, "path"), repfile), 
                                         row.names = 1)))
  isMature = size >= .getPar(this, "species.maturity.size")
  
  MLF = .calculateMLF(fecundity, isMature, weight)
  
  return(MLF)
  
}


calculateMortality = function(conf, sp) {
  
  .getM = function(n, ratio, MLF, tn, d1, value=TRUE) {
    G = -log((1/ratio)/MLF)/n
    alpha = (tn/d1)^(1/(n-1)) - 1
    di = c(d1, alpha*(1+alpha)^(2:n - 2)*d1)
    Mi = G/di
    if(isTRUE(value)) return(list(Mi=Mi, di=di))
    return(all(diff(Mi) < 0))
  }
  
  this = .getPar(conf, sp=sp)
  d1 = .getPar(this, "species.egg.stage.duration") # days
  if(is.null(d1)) d1 = 2
  d1 = d1/365 # transformed to years
  
  tsMLF = calculateMLF(conf, sp=sp)
  MLF = mean(tsMLF) 
  
  this = .getPar(conf, sp=sp)
  
  tn = .getPar(this, "species.lifespan")
  ratio = .getPar(this, "species.sexratio")
  
  ind = sapply(2:20, .getM, ratio=ratio, MLF=MLF, tn=tn, d1=d1, value=FALSE)
  n = which.min(ind)
  
  tmp = .getM(n=n, ratio=ratio, MLF=MLF, tn=tn, d1=d1)
  Mi = tmp$Mi
  di = tmp$di
  
  ti = cumsum(di)
  
  size = VB(c(0, ti), this, method=3)
  size[length(size)] = Inf
  out = list(age = c(0, ti), size=size, M=Mi)
  return(out)
  
}





# Internal ----------------------------------------------------------------

read.cal = function(conf, sp) {
  
  this = .getPar(conf, sp=sp)
  ndt = conf$simulation.time.ndtperyear
  T   = conf$simulation.time.nyear*ndt
  
  landings = read.yield(conf, sp)
 
  nonHarvest = all(landings==0)
  
  if(isTRUE(nonHarvest)) {
    
    Linf = .getPar(this, "species.linf")
    bins = pretty(c(0, 0.9*Linf), n=15)
    dbin = unique(diff(bins))
    length_classes = 0.5*head(bins, -1) + 0.5*tail(bins, -1)
    
    newmat = matrix(0, nrow=T, ncol=length(length_classes))
    
    output = list(cal=NULL, marks=length_classes, dbin=dbin, mat=newmat, bins=bins,
                  nonHarvest=TRUE)
    return(output)
    
  }
  
  a = .getPar(this, "species.length2weight.condition.factor")
  b = .getPar(this, "species.length2weight.allometric.power")
  
  file = .getPar(this, "fisheries.catchatlength.file")
  if(length(file)>1) stop("Only one catch-at-length file must be provided.")
  if(is.null(file)) return(NULL)
  
  file = file.path(attr(file, "path"), file)
  periods = c("year", "quarter", "month", "week")
  out = read.csv(file, check.names = FALSE)
  must = names(out)[names(out) %in% periods]
  if(length(must)<1) stop("Missing time information in catch-at-length file.")
  check = !(must %in% names(out))
  if(any(check)) stop("Missing time information in catch-at-length file.")
  length_classes = as.numeric(setdiff(colnames(out), must))
  bad = paste(setdiff(colnames(out), must)[is.na(length_classes)], collapse=", ")
  if(any(is.na(length_classes)))
    stop(sprintf("Size class marks should be numeric, check: %s.", bad))
  check = !identical(length_classes, sort(length_classes))
  if(check) stop("Catch-at-length size classes must be in increasing order.")
  dbin = unique(diff(length_classes))
  mat = as.matrix(out[, as.character(length_classes)])

  ndtcal = .getPar(this, "fisheries.catchatlength.ndtPerYear")
  if(is.null(ndtcal)) stop("Parameter 'fisheries.catchatlength.ndtPerYear' is missing.")  
  
  ix = .time.conv(ndtcal, ndt, nrow(mat), T)
  bins = c(length_classes - dbin, length_classes[length(length_classes)] + dbin)

  isize = pmax(bins, .getPar(this, "egg.size")) 
  
  L1 = head(isize, -1)
  L2 = tail(isize, -1)
  W2 = (a/(b+1))*(L2^(b+1) - L1^(b+1))
  
  newmat = ix$w*mat[ix$ind, ]
  
  wmat = t(t(newmat)*W2)
  ilandings = 1e-6*rowSums(wmat)
  ilandings[ilandings==0] = 1
  units = landings/ilandings
  if(all(is.na(units))) {
    units = 1 # assume CAL is unbiased
    warning("No landing data, assuming catch-at-length is unbiased and representing the full landings. 
            Manually calculate the landings in recommended.")
  }
  
  units[is.na(units)] = 1 # assume unbiased when landing data is not available.
  
  newmat = newmat*units
  
  output = list(cal=out, marks=length_classes, dbin=dbin, mat=newmat, bins=bins, 
                nonHarvest=FALSE)
  return(output)
}

rebinning = function(x, y) {
  
  if(is.list(y)) stop("'y' must be a single vector.")
  
  .rebinning = function(x, y, k=10000) {
    .mini = function(x, k=100) head(approx(x=x, n=k*length(x)-(k-1))$y, -1)
    .mytable = function(x, levels) table(factor(x, levels=levels))
    xm = cbind(head(x, -1), tail(x, -1))
    mini = apply(xm, 1, .mini, k=k)
    out = cut(mini, breaks=y, right = FALSE, include.lowest=TRUE)
    levels = levels(out)
    out = matrix(out, ncol=k, byrow = TRUE)
    xout = t(apply(out, 1, .mytable, levels=levels))/k
    rownames(xout) = levels(cut(x, breaks = x, right = FALSE, include.lowest=TRUE))
    return(xout)
  }
  
  if(!is.list(x)) return(.rebinning(x, y))
  
  return(lapply(x, .rebinning, y=y))
  
}


.time.conv = function(ndt_in, ndt, Tref, T) {
  caltime = seq(from=0, by=1/ndt_in, length.out=Tref+1)
  simtime = seq(from=0.5/ndt, by=1/ndt, length.out=T)
  ind = cut(simtime, breaks = caltime, labels = FALSE)
  const = rle(ind)
  const$values = 1/const$lengths
  w = inverse.rle(const)
  return(list(ind=ind, w=w))
}

read.biomass = function(conf, sp) {
  
  this = .getPar(conf, sp=sp)
  ndt = conf$simulation.time.ndtperyear 
  T = ndt*conf$simulation.time.nyear
  biofile = .getPar(this, "observed.biomass.file")
  if(is.null(biofile)) stop("Observed biomass have not been provided.")
  bioref = .readCSV(biofile)
  ivar= .getPar(this, "species.name")
  ndtbio = .getPar(this, "observed.biomass.ndtPerYear")
  if(is.null(ndtbio)) stop("Parameter 'observed.biomass.ndtPerYear' is missing.")
  ix = .time.conv(ndtbio, ndt, nrow(bioref), T)
  biomass = bioref[ix$ind, ivar]
  return(biomass)
  
}

read.yield = function(conf, sp) {
  
  this = .getPar(conf, sp=sp)
  ndt = conf$simulation.time.ndtperyear 
  T = ndt*conf$simulation.time.nyear
  biofile = .getPar(this, "fisheries.yield.file")
  if(is.null(biofile)) stop("Landings have not been provided.")
  bioref = .readCSV(biofile)
  ivar= .getPar(this, "species.name")
  ndtbio = .getPar(this, "fisheries.yield.ndtPerYear")
  if(is.null(ndtbio)) stop("Parameter 'fisheries.yield.ndtPerYear' is missing.")
  ix = .time.conv(ndtbio, ndt, nrow(bioref), T)
  biomass = ix$w*bioref[ix$ind, ivar]
  return(biomass)
  
}

read.fecundity = function(conf, sp) {
  
  this = .getPar(conf, sp=sp)
  repfile = .getPar(this, "reproduction.season.file")
  fecundity = as.numeric(unlist(.readCSV(repfile, row.names = 1)))
  return(fecundity)
  
}


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
  nonHarvest = CAL$nonHarvest
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
                  bins=list(age=age_bins, size=size_bins), nonHarvest=TRUE)
    
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
    
    biopop = 1e-6*t(t(pop)*weight)[, -c(1:xn)]
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
  pind = xind + 0:23
  
  es = colSums(catch[pind, ])/colMeans(pop[pind, ])
  
  Fguess = sum(yield[pind]/biomass[pind])
  
  output = list(pop=pop[pind, ], catch=catch[pind, ], R=R, 
                biomass=biomass[pind], yield=yield[pind], F=Fseason,
                dist = pop[xind ,], distB = pop[xind ,]*weight, 
                selectivity=es, age=age, size=size,
                Fguess=Fguess, observed=list(biomass=bioguess, yield=yield_obs[1:ndt]),
                bins=list(age=age_bins, size=size_bins), nonHarvest=nonHarvest)
  
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
  nonHarvest = if(is.null(CAL)) FALSE else CAL$nonHarvest
  if(isTRUE(nonHarvest)) return(.simCatch_ini(conf, sp))
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
    
    biopop = 1e-6*t(t(pop)*weight)[, -c(1:xn)]
    biomass = rowSums(biopop)
    yield   = 1e-6*rowSums(t(t(catch)*weight))

    ll_biomass  = lnorm2(biofit, tail(biomass, ndt))
    ll_yield    = lnorm2(yield_obs[1:ndt] , tail(yield, ndt))
    
    ll = ll_biomass*llw(cv[1]) + ll_yield*llw(cv[2])
    
    if(!isTRUE(value)) return(ll)

    xind = T - ndt + 1
    pind = xind + 0:23
    
    es = colSums(catch[pind, ])/colMeans(pop[pind, ])
    
    output = list(pop=pop[pind, ], catch=catch[pind, ], R=R, 
                  biomass=biomass[pind], yield=yield[pind], F=Fseason, Fguess=F,
                  dist=pop[xind ,], distB = 1e-6*pop[xind ,]*weight, 
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

.bioguess = function(x, ndt, ts=FALSE) {
  if(all(is.na(x))) stop("No biomass information is provided.") # check
  out_ts = x[1:ndt]
  out = mean(out_ts, na.rm=TRUE)
  if(is.na(out)) {
    warning("No biomass information for the first year, using average of time series.")
    out = mean(x, na.rm=TRUE)
    out_ts = rep(NA, ndt)
    out_ts[1] = out
  }
  if(isTRUE(ts)) return(out_ts) else return(out)
} 

.getSelectivity = function(size, this) {
  
  par = list(type = .getPar(this, "fisheries.selectivity.type"),
             L50  = .getPar(this, "fisheries.selectivity.l50"),
             L75  = .getPar(this, "fisheries.selectivity.l75"),
             tiny = .getPar(this, "fisheries.selectivity.tiny"))
  
  if(is.null(par$tiny)) par$tiny = 1e-3
  
  return(.calculateSelectivity(x=size, par=par))
  
}

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

.readCSV = function(file, ...) read.csv(file.path(attr(file, "path"), file), ...)


# OSMOSE package ----------------------------------------------------------


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
    sim = .simF_ini(conf, sp)
    sim$osmose = .initial_length_dist(sim, sp)
    pars = rbind(pars, as.matrix(sim$osmose))
    this = .getPar(conf, sp=sp)
    out[[.getPar(this, "species.name")]] = sim
  }
  
  pars = as.data.frame(pars)
  colnames(pars) = NULL
  
  xoutput = list(par=pars, init=out)
  class(xoutput) = c("osmose.initialization", class(xoutput))
  
  if(!is.null(output)) osmose::write_osmose(pars, file=output, sep=" = ")
    
  return(invisible(xoutput))
  
}

