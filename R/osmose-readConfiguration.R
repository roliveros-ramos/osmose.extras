

# Main --------------------------------------------------------------------


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


# Internal ----------------------------------------------------------------

.getPar = function(conf, par, sp=NULL) {
  if(!is.null(sp)) par = sprintf(".sp%d$", sp)
  par = tolower(par)
  out = conf[grep(names(conf), pattern=par)]
  if(length(out)==1) out = out[[1]]
  if(length(out)==0) out = NULL
  return(out)
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

.readCSV = function(file, ...) {
  
  .guessSeparator = function(Line) {
    SEPARATORS = c(equal = "=", semicolon = ";",
                   coma = ",", colon = ":", tab = "\t")
    guess = which.min(nchar(lapply(str_split(Line,SEPARATORS), "[", i = 1)))
    separator = SEPARATORS[guess]
    
    return(separator)
  }
  
  config = readLines(file) # read lines
  sep = names(which.max(table(sapply(config, .guessSeparator))))
  
  file = file.path(attr(file, "path"), file)
  
  out = if(sep==",") read.csv(file=file, ...) else read.csv2(file=file, sep=sep, ...)
  return(out)
  
}

# Data Parsing ------------------------------------------------------------

.setupInitialization = function(conf) {
  
  nsp = .getPar(conf, "simulation.nspecies")
  
  spind = .getPar(conf, "species.type") == "focal"
  spind = gsub(names(spind)[which(spind)], pattern="species.type.sp", replacement = "") 
  spnames = .getPar(conf, "species.name")[sprintf("species.name.sp%s", spind)]
  spind = as.numeric(spind)
  
  for(sp in spind) {
    
    sim = list()
    sim$cal       = read.cal(conf, sp)
    sim$biomass   = read.biomass(conf, sp)
    sim$yield     = read.yield(conf, sp)
    sim$fecundity = read.fecundity(conf, sp)
    isp = sprintf("osmose.initialization.data.sp%d", sp)
    conf[[isp]]   = sim
    
  }
  
  return(conf)
  
}


read.cal = function(conf, sp) {
  
  this = .getPar(conf, sp=sp)
  ndt = conf$simulation.time.ndtperyear
  T   = conf$simulation.time.nyear*ndt
  
  spname = .getPar(this, "species.name")
  landings = read.yield(conf, sp)
  
  harvested = !all(landings==0)
  
  if(!isTRUE(harvested)) {
    
    Linf = .getPar(this, "species.linf")
    bins = pretty(c(0, 0.9*Linf), n=15)
    dbin = unique(diff(bins))
    length_classes = 0.5*head(bins, -1) + 0.5*tail(bins, -1)
    
    newmat = matrix(0, nrow=T, ncol=length(length_classes))
    
    output = list(cal=NULL, marks=length_classes, dbin=dbin, mat=newmat, bins=bins,
                  harvested=FALSE)
    return(output)
    
  }
  
  a = .getPar(this, "species.length2weight.condition.factor")
  b = .getPar(this, "species.length2weight.allometric.power")
  
  file = .getPar(this, "fisheries.catchatlength.file")
  msg = sprintf("Only one catch-at-length file must be provided for %s.", spname)
  if(length(file)>1) stop(msg)
  if(is.null(file)) return(NULL)
  
  file = file.path(attr(file, "path"), file)
  periods = c("year", "quarter", "month", "week")
  out = read.csv(file, check.names = FALSE)
  must = names(out)[names(out) %in% periods]
  msg = sprintf("Missing time information in %s's catch-at-length file.", spname)
  if(length(must)<1) stop(msg)
  check = !(must %in% names(out))
  if(any(check)) stop(msg)
  length_classes = as.numeric(setdiff(colnames(out), must))
  bad = paste(setdiff(colnames(out), must)[is.na(length_classes)], collapse=", ")
  if(any(is.na(length_classes)))
    stop(sprintf("Size class marks should be numeric, check: %s in %s's file.", bad, spname))
  check = !identical(length_classes, sort(length_classes))
  msg = sprintf("Catch-at-length size classes must be in increasing order, check %s's file.", spname)
  if(check) stop(msg)
  dbin = unique(diff(length_classes))
  
  mat = as.matrix(out[, as.character(length_classes)])

  ndtcal = .getPar(this, "fisheries.catchatlength.ndtPerYear")
  msg = sprintf("Parameter 'fisheries.catchatlength.ndtPerYear.sp%d' is missing.", sp)
  if(is.null(ndtcal)) stop(msg)  
  nT = conf$simulation.time.nyear*ndtcal 
  
  msg = sprintf("Catch-at-length data for %s is incomplete, %d rows found and at least %d rows are expected (%d x %d).",
                spname, nrow(mat), nT, conf$simulation.time.nyear, ndtcal)
  if(nrow(mat)<nT) stop(msg)
  msg = sprintf("Only first %d rows of %s's catch-at-length data are used.", nT, spname)
  if(nrow(mat)>nT) warning(msg)
  mat = mat[seq_len(nT), ]
  
  if(all(is.na(mat))) {
    
    msg = sprintf("No catch-at-length data in %s's file is provided, all NAs.", spname)
    warning(msg)
    
    Linf = .getPar(this, "species.linf")
    bins = pretty(c(0, 0.9*Linf), n=15)
    dbin = unique(diff(bins))
    length_classes = 0.5*head(bins, -1) + 0.5*tail(bins, -1)
    # check for the right filling
    newmat = matrix(NA, nrow=T, ncol=length(length_classes))
    output = list(cal=NULL, marks=length_classes, dbin=dbin, mat=newmat, bins=NULL,
                  harvested=TRUE)
    
  }
  
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
    msg = sprintf("No landing data for %s, assuming catch-at-length is unbiased and representing the full landings. 
            If TRUE, please manually calculate the landings from catch-at-length.", spname)
    warning(msg)
  }
  
  units[is.na(units)] = 1 # assume unbiased when landing data is not available.
  
  newmat = newmat*units
  
  output = list(cal=out, marks=length_classes, dbin=dbin, mat=newmat, bins=bins, 
                harvested=TRUE)
  return(output)
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
  sfec = sum(fecundity) 
  if(abs(sfec-1)< 1e-2) {
    warning("Using relative fecundities.")
    relfec = .getPar(this, "species.relativefecundity")
    fecundity = relfec*fecundity
  }
  return(fecundity)
  
}


