
####################
### PLM EXTENSION - STATA-LIKE DEGREE OF FREEDOM ADJUSTMENT FOR DIFFERENT TYPES OF
### CORRECTED STANDARD ERRORS
####################

# for within transformation ONLY
# same df correction applies in case of 2SLS (i.e k equals thr number of model regressors NOT 
# the number of exogenous regressors (incl. instruments))

######################################################

### HETEROSKEDASTICITY CONSISTENT (WHITE'S) STANDARD ERRORS

vcovHR = function(x, stata = T) {
  
  demy = pmodel.response(x, model = "within")
  demX = model.matrix(x, model = "within", rhs = 1, cstcovar.rm = "all")
  
  if (length(formula(x))[2] > 1) { # check for IV option
    demZ = model.matrix(x, model = "within", rhs = 2, cstcovar.rm = "all")
    demX = fitted(lm.fit(demZ, demX))%>%as.matrix()
  }
  
  pdim = pdim(x)
  nT = pdim$nT$N # nmb of obs
  k = dim(demX)[[2]] # nmb of x-vars
  n = pdim$nT$n # nmb of groups
  t = pdim$nT$T # nmb of time units
  
  k_ef = case_when(
    x$args$effect  == "individual" ~ as.integer(n),
    x$args$effect  == "time" ~ as.integer(t),
    x$args$effect  == "twoways" ~ as.integer((n+t-1)) ) # cannot substract intercept twice
  
  uhat = x$residuals
  
  salame = crossprod(demX, diag(uhat^2)) %*% demX 
  pane = solve(crossprod(demX))
  mycov = pane %*% salame %*% pane
  
  if (stata) {
    mycov = nT/(nT-k-k_ef) * mycov
  }
  
  return(mycov)
  
}

######################################################

### ONEWAY CLUSTERED (ROGERS') STANDARD ERRORS

# higher-level clustering possible
# remove missing values FIRST
# order by group-time FIRST

vcovCL = function(x, cluster, stata = T) {
  
  demy = pmodel.response(x, model = "within")
  demX = model.matrix(x, model = "within", rhs = 1, cstcovar.rm = "all")
  
  if (length(formula(x))[2] > 1) { # check for IV option
    demZ = model.matrix(x, model = "within", rhs = 2, cstcovar.rm = "all")
    demX = fitted(lm.fit(demZ, demX))%>%as.matrix()
  }
  
  pdim = pdim(x)
  nT = pdim$nT$N # nmb of obs
  k = dim(demX)[[2]] # nmb of x-vars
  n = pdim$nT$n # nmb of groups
  t = pdim$nT$T # nmb of time units
  glist = cluster %>% unique()
  g = glist %>% length() # nmb of cluster groups
  
  k_ef = case_when(
    x$args$effect  == "individual" ~ as.integer(n),
    x$args$effect  == "time" ~ as.integer(t),
    x$args$effect  == "twoways" ~ as.integer((n+t-1)) ) # cannot substract intercept twice
  
  uhat = x$residuals
  
  selector = vector("list", g)
  for (i in 1:g) {
    selector[[i]] = which(cluster == glist[i])
  }
  
  Sl = array(dim = c(k, k, g))
  for (i in 1:g) {
    X = demX[selector[[i]], , drop = FALSE]
    u = uhat[selector[[i]]]
    Sl[, , i] = crossprod(X, outer(u, u)) %*% X
  }
  salame = apply(Sl, 1:2, sum)
  pane = solve(crossprod(demX))
  mycov = pane %*% salame %*% pane
  
  if (stata) {
    mycov = (nT-1)/(nT-k-k_ef) * g/(g-1) * mycov
  }
  
  return(mycov)
  
}

######################################################

### TWOWAY CLUSTERED (ROGERS') STANDARD ERRORS

# higher-level clustering possible
# remove missing values FIRST
# order by group-time FIRST

vcovTC = function(x, cluster1, cluster2, stata = T) {
  
  mycov = vcovCL(x, cluster = cluster1, stata = stata) + 
    vcovCL(x, cluster = cluster2, stata = stata) -
    vcovHR(x, stata = stata)
  
  return(mycov)
    
}
