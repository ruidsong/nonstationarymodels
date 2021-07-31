glmssn_nonstationary <- function(formula, 
                                 ssn.object, 
                                 dichvar,
                                 family = "Gaussian",
                                 CorModels = NULL,
                                 use.nugget = TRUE,
                                 use.anisotropy = FALSE,
                                 addfunccol = NULL,
                                 trialscol = NULL,
                                 EstMeth = "REML",
                                 useTailDownWeight = FALSE,
                                 trans.power = NULL,
                                 trans.shift = 0,
                                 control = list(max.range.factor = 4,
                                                trunc.pseudo = NULL,
                                                maxiter.pseudo = 20,
                                                beta.converge = 1e-5), 
                                 init,
                                 upper.bound = c(Inf, Inf, Inf, Inf, Inf, Inf),
                                 lower.bound = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf)){ 
  Warnlog <- NULL
  data <- ssn.object@obspoints@SSNPoints[[1]]@point.data
  data <- cbind(data, ssn.object@obspoints@SSNPoints[[1]]@point.coords)
  xcol <- "coords.x1"
  ycol <- "coords.x2"
  family <- tolower(family)
  
  if(!"max.range.factor" %in% names(control)) control$max.range.factor <- 4
  if(!"maxiter.pseudo" %in% names(control)) control$maxiter.pseudo <- 20
  
  nIDs <- sort(as.integer(as.character(unique(data[,"netID"]))))
  dist.junc <- matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
  net.zero <-  matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
  nsofar <- 0
  distord <- order(data[,"netID"],data[,"pid"])
  names(distord) <- rownames(data)[distord]
  for(i in nIDs) {
    workspace.name <- paste("dist.net", i, ".RData", sep = "")
    path <- file.path(ssn.object@path, "distance", "obs",
                      workspace.name)
    if(!file.exists(path)) {
      stop("Unable to locate required distance matrix")
    }
    file_handle <- file(path, open="rb")
    distmat <- unserialize(file_handle)
    ordpi <- order(as.numeric(rownames(distmat)))
    close(file_handle)
    ni <- length(distmat[1,])
    dist.junc[(nsofar + 1):(nsofar + ni),(nsofar + 1):
                (nsofar + ni)] <- distmat[ordpi, ordpi, drop = F]
    net.zero[(nsofar + 1):(nsofar + ni),(nsofar + 1):(nsofar + ni)] <- 1
    nsofar <- nsofar + ni
  }
  
  mf <- match.call(expand.dots = FALSE)
  dataXY.out <- dataXY(formula, data,
                       family = family, trialscol = trialscol,
                       trans.power = trans.power,
                       trans.shift = trans.shift,
                       CorModels = CorModels,
                       distord = distord)
  
  REs <- dataXY.out$REs
  REmodelmatrices <- dataXY.out$REmodelmatrices
  n.all <- dataXY.out$sampsizes$n.all
  a.mat <- NULL
  b.mat <- NULL
  a.mat.data <- NULL
  b.mat.data <- NULL
  dist.hydro <- NULL
  dist.hydro.data <- NULL
  w.matrix.data <- NULL
  ind <- dataXY.out$indvecs$ind.allxy
  
  ## create any necessary matrices from distance and flow matrices
  
  if(!is.null(dist.junc) ) {
    ## maximum distance to common junction between two sites
    a.mat <- pmax(dist.junc,t(dist.junc))
    a.mat.data <- a.mat[ind,ind]
    ## minimum distance to common junction between two sites
    b.mat <- pmin(dist.junc,t(dist.junc))
    b.mat.data <- b.mat[ind,ind]
    ## hydrological distance
    dist.hydro <- as.matrix(dist.junc + t(dist.junc))
    ## subset stream distance to observed locations only
    dist.hydro.data <- dist.hydro[ind, ind]
  }
  
  if(length(grep("tailup",CorModels)) | useTailDownWeight == TRUE) {
    if(missing(addfunccol) || is.null(addfunccol) ||
       length(addfunccol) == 0 || !(addfunccol %in% colnames(data)))
      stop("The specified value for addfunccol was invalid")
    flow.con.mat <- 1 - (b.mat > 0)*1
    
    w.matrix <- sqrt(pmin(outer(data[distord,addfunccol], rep(1, times = n.all)),
                          t(outer(data[distord,addfunccol],rep(1, times = n.all)))) /
                       pmax(outer(data[distord,addfunccol],rep(1, times = n.all)),
                            t(outer(data[distord,addfunccol],rep(1, times = n.all)))))*
      flow.con.mat * net.zero
    w.matrix.data <- w.matrix[ind, ind]
  }
  net.zero.data <- net.zero[ind,ind]
  
  xcoord <- data[distord,xcol]
  ycoord <- data[distord,ycol]
  xcoord.data <- xcoord[ind]
  ycoord.data <- ycoord[ind]
  
  z <- dataXY.out$respvecs$z
  X2 <- dataXY.out$Xmats$X2
  n.allxy <- dataXY.out$sampsizes$n.allxy
  
  trialsvec <- NULL
  ## Initial parameter estimates, fixed effects, pseudo-data
  if(family == "gaussian"){
    A.5 <-  NULL
    Del.i <-  NULL
    zt <- z
  }
  
  ## We're going to pass this environment in to the m2LL.stream
  ## funtion and it will track all evaluations of the -loglik surface
  loglik.environment <- environment()
  assign("RESULT",NULL,loglik.environment)
  
  # set maximum range parameters
  maxrang <- NULL
  mrf <- control$max.range.factor
  if(length(grep("tailup",CorModels)) > 0)
    maxrang <- c(maxrang, NA, mrf*max(dist.hydro))
  if(length(grep("taildown",CorModels)) > 0)
    maxrang <- c(maxrang, NA, mrf*max(dist.hydro))
  if(length(grep("Euclid",CorModels)) > 0)
    maxrang <- c(maxrang, NA, mrf*max(distGeo(xcoord.data,
                                              ycoord.data, xcoord.data, ycoord.data, 1)))
  if(length(REs)) maxrang <- c(maxrang, rep(NA, times = length(REs)))
  if(use.nugget == TRUE) maxrang <- c(maxrang, NA)
  
  # ----------- START LOOPING HERE ---------------------------
  
  # create an indicator to stop looping
  stoploop <- 0
  # keep track of the number of iterations
  iter <- 0
  # keep track of number of inner iterations for beta
  inner.iter2 <- NULL
  # begin looping
  while(stoploop == 0) {
    ## initial parameter estimates
    if(iter == 0) {
      theta <- theta.ini(z = zt, X = X2,
                         CorModels= CorModels,
                         use.nugget = use.nugget, use.anisotropy = use.anisotropy,
                         dist.hydro.data = dist.hydro.data, x.dat = xcoord.data,
                         y.dat = ycoord.data, REs = REs, init = init)
      attributes(theta) ## scale, type, terms
      TH.scale <- attributes(theta)$scale
      TH.type <- attributes(theta)$type
      TH.terms <- attributes(theta)$terms
    }
    
    if(length(theta > 0)) {
      ## optimizing for covariance parameter estimates using ML or REML
      if(length(theta) == 1) {
        lowerb <- log(.1*exp(theta))
        upperb <- log(10*exp(theta))
        parmest1.out <- optimize(m2LL.stream, interval = c(lowerb,upperb), 
                                 m2LLdata = zt, X = X2,
                                 dist.hydro = dist.hydro.data, weight = w.matrix.data,
                                 net.zero = net.zero.data,
                                 a.mat = a.mat.data, b.mat = b.mat.data,
                                 x.dat = xcoord.data, y.dat = ycoord.data,
                                 Del.i = Del.i, A.5 = A.5,
                                 CorModels = CorModels, scale = TH.scale, 
                                 useTailDownWeight = useTailDownWeight,
                                 use.nugget = use.nugget, use.anisotropy = use.anisotropy,
                                 EstMeth = EstMeth, loglik.environment=loglik.environment,
                                 REs = REs, maxrang = maxrang)
        lowerb <- log(.5*exp(parmest1.out$minimum))
        upperb <- log(2*exp(parmest1.out$minimum))
        parmest2.out <- optimize(m2LL.stream, interval =
                                   c(lowerb,upperb), m2LLdata = zt, X = X2,
                                 dist.hydro = dist.hydro.data, weight = w.matrix.data,
                                 net.zero = net.zero.data,
                                 a.mat = a.mat.data, b.mat = b.mat.data,
                                 x.dat = xcoord.data, y.dat = ycoord.data,
                                 Del.i = Del.i, A.5 = A.5,
                                 CorModels = CorModels, scale = TH.scale, 
                                 useTailDownWeight = useTailDownWeight,
                                 use.nugget = use.nugget, use.anisotropy = use.anisotropy,
                                 EstMeth = EstMeth, loglik.environment=loglik.environment,
                                 REs = REs, maxrang = maxrang)
        parmest.out <- parmest2.out
        theta <- parmest2.out$minimum
        m2LL <- parmest2.out$objective
      }
      
      if(length(theta) > 1) {
        if(iter == 0) {
          parmest1.out <- optim(theta, m2LL.stream, m2LLdata = zt, dichvar = dichvar,
                                X = X2, dist.hydro = dist.hydro.data,
                                weight = w.matrix.data, net.zero = net.zero.data,
                                a.mat = a.mat.data, b.mat = b.mat.data,
                                x.dat = xcoord.data, y.dat = ycoord.data,
                                Del.i = Del.i, A.5 = A.5,
                                CorModels = CorModels, useTailDownWeight = useTailDownWeight,
                                use.nugget = use.nugget, use.anisotropy = use.anisotropy,
                                EstMeth = EstMeth, method = "L-BFGS-B", upper = upper.bound, lower = lower.bound,
                                loglik.environment=loglik.environment, REs = REs, 
                                scale = TH.scale, maxrang = maxrang, data = data, control = list(maxit = 1000, factr = 1e-5))
        }
        
        parmest2.out <- optim(theta, m2LL.stream, m2LLdata = zt, dichvar = dichvar,
                              X = X2, dist.hydro = dist.hydro.data, weight = w.matrix.data,
                              net.zero = net.zero.data,
                              a.mat = a.mat.data, b.mat = b.mat.data,
                              x.dat = xcoord.data, y.dat = ycoord.data,
                              Del.i = Del.i, A.5 = A.5,
                              CorModels = CorModels, useTailDownWeight = useTailDownWeight,
                              use.nugget = use.nugget, use.anisotropy = use.anisotropy,
                              EstMeth = EstMeth, method = "L-BFGS-B", upper = upper.bound, lower = lower.bound,
                              loglik.environment=loglik.environment,
                              scale = TH.scale, REs = REs, maxrang = maxrang, data = data, control = list(maxit = 1000, factr = 1e-5))
        if(iter == 0 & parmest1.out$value < parmest2.out$value) {
          parmest.out <- parmest1.out
        } else parmest.out <- parmest2.out
        theta <- parmest.out$par
        m2LL <- parmest.out$value
      }
      
      ## go back to original scale for covariance parameters
      parmest <- untrans.theta(theta = theta, scale = TH.scale)
    }
    
    if(is.null(theta)) parmest <- 1
    if(length(CorModels) == 0) V <- diag(rep(parmest, times = n.allxy))
    else V <- makeCovMat_nonstationary(parmest, dichvar = dichvar, 
                         dist.hydro = dist.hydro.data,
                         w.matrix = w.matrix.data, net.zero = net.zero.data,
                         a.mat = a.mat.data, b.mat = b.mat.data, 
                         x.row = xcoord.data, y.row = ycoord.data,
                         x.col = xcoord.data, y.col = ycoord.data,
                         CorModels = CorModels, useTailDownWeight = useTailDownWeight,
                         use.nugget = use.nugget, use.anisotropy = use.anisotropy, REs = REs, data = data)
    
    if(!is.null(Del.i)) {
      ## V <- diag(Del.i) %*% diag(A.5) %*% V %*% diag(A.5) %*% diag(Del.i)
      V <- Del.i*A.5*t((Del.i*A.5) * V)
    }
    ##	if(min(svd(V)$d) < 0 ) {browser()}
    qrV <- qr(V)
    ##	if(class(qrV) != "qr") {browser()}
    ViX <- solve(qrV,X2)
    covbi <- t(X2) %*% ViX
    covb <- solve(covbi)
    if(family != "gaussian") beta.old <- beta.hat
    beta.hat <- covb %*% t(ViX) %*% zt
    if(family == "gaussian") beta.old <- beta.hat
    if(family == "gaussian"){
      A.5 <-  NULL
      Del.i <-  NULL
      zt <- z
    }
    
    ## convergence criteria on the fixed effect parameters
    non0ind <- beta.old != 0
    if(all(abs((beta.hat[non0ind] - beta.old[non0ind])/
               beta.old[non0ind]) < control$beta.converge)) stoploop <- 1
    if (iter > control$maxiter.pseudo) {
      Warnlog <- c(Warnlog,
                   paste("More than", control$maxiter.pseudo,
                         "iterations, algorithm terminated -- results questionable"))
      stoploop <- 1
    }
    iter <- iter + 1
  }
  # ----------- DONE LOOPING HERE ---------------------------
  
  # inverse covariance matrix between observed locations
  Vi <- solve(qrV)
  ViX <- solve(qrV,X2)
  covbi <- t(X2) %*% ViX
  covb <- solve(covbi)
  bhat.se <- sqrt(diag(covb))
  b.hat <- covb %*% t(X2) %*% Vi %*% zt
  bhat.se <- sqrt(diag(covb))
  p <-  dataXY.out$Xmats$p
  if(use.nugget == TRUE)  nugget <- parmest[length(parmest)]
  is.na(data[,dataXY.out$respvecs$response.col])
  nobs <- length(ssn.object@obspoints@SSNPoints[[1]]@point.data[,1])
  
  # if any missing data, create prediction set in the SSN object
  if(dataXY.out$sampsizes$n.allcov > n.allxy) {
    TempSSNPoints <- ssn.object@obspoints@SSNPoints
    TempPredPoints <- ssn.object@predpoints
    listlen <- length(TempPredPoints@SSNPoints)
    TempPredPoints@SSNPoints[[listlen+1]] <- TempSSNPoints[[1]]
    TempPredPoints@SSNPoints[[listlen+1]]@network.point.coords <-
      TempPredPoints@SSNPoints[[listlen+1]]@network.point.coords[
        is.na(data[,dataXY.out$respvecs$response.col]), ,drop = F]
    TempPredPoints@SSNPoints[[listlen+1]]@point.coords <-
      junk <- TempPredPoints@SSNPoints[[listlen+1]]@point.coords[
        is.na(data[,dataXY.out$respvecs$response.col]), ,drop = F]
    TempPredPoints@SSNPoints[[listlen+1]]@point.data <-
      TempPredPoints@SSNPoints[[listlen+1]]@point.data[
        is.na(data[,dataXY.out$respvecs$response.col]), ,drop = F]
    TempPredPoints@ID <- c(TempPredPoints@ID, "_MissingObs_")
    ssn.object@predpoints <- TempPredPoints
  }
  
  if(is.null(theta)) {
    nugget <- 1
    parmest <- 1
    attr(parmest,"scale") <- "natural"
    attr(parmest, "type") <- "parsill"
    attr(parmest, "terms") <- "nugget(held at 1)"
    m2LL <- NULL
    parmest.out <- NULL
  } else {
    attr(parmest,"scale") <- TH.scale
    attr(parmest,"type") <- TH.type
    attr(parmest,"terms") <- TH.terms
  }
  
  if(use.nugget == FALSE) nugget <- 0
  outpt <- list(
    args = list(
      formula = formula,
      zcol = dataXY.out$respvecs$response.col,
      family = family,
      CorModels = CorModels,
      useTailDownWeight = useTailDownWeight,
      use.nugget = use.nugget,
      use.anisotropy = use.anisotropy,
      addfunccol = addfunccol,
      trialscol = trialscol,
      EstMeth = EstMeth,
      trans.power = trans.power,
      trans.shift = trans.shift,
      algorithm = "orig"),
    
    ssn.object = ssn.object,
    
    sampinfo = list(
      ind.obs = ind,
      ind.RespNA = dataXY.out$indvecs$ind.RespNA,
      sample.size = nobs,
      obs.sample.size = n.allxy,
      missing.sample.size = nobs - n.allxy,
      rankX = p,
      z = zt,
      trialsvec = trialsvec,
      X = X2,
      effnames = dataXY.out$Xmats$effnames,
      setzero = dataXY.out$indvecs$setzero,
      setNA = dataXY.out$indvecs$setNA,
      setNA2 = dataXY.out$indvecs$setNA2,
      cutX1toX2 = dataXY.out$indvecs$cutX1toX2,
      REs = REs,
      REmodelmatrices = REmodelmatrices),
    
    estimates = list(
      theta = parmest,
      nugget = nugget,
      V = V,
      Vi = Vi,
      betahat = b.hat,
      covb = covb,
      covbi = covbi,
      m2LL = m2LL,
      Warnlog = Warnlog),
    
    loglik.surface=get("RESULT",loglik.environment),
    
    optimOutput=parmest.out
  )
  
  class(outpt) <- "glmssn"
  outpt
}



grid_search <- function(formula, 
                        ssn.object, 
                        dichvar,
                        family = "Gaussian",
                        CorModels = NULL,
                        use.nugget = TRUE,
                        use.anisotropy = FALSE,
                        addfunccol = NULL,
                        trialscol = NULL,
                        EstMeth = "REML",
                        useTailDownWeight = FALSE,
                        trans.power = NULL,
                        trans.shift = 0,
                        upper.bound,
                        lower.bound,
                        grid_depth) {
  data <- ssn.object@obspoints@SSNPoints[[1]]@point.data
  data <- cbind(data, ssn.object@obspoints@SSNPoints[[1]]@point.coords)
  xcol <- "coords.x1"
  ycol <- "coords.x2"
  family <- tolower(family)
  
  nIDs <- sort(as.integer(as.character(unique(data[,"netID"]))))
  dist.junc <- matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
  net.zero <-  matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
  nsofar <- 0
  distord <- order(data[,"netID"],data[,"pid"])
  names(distord) <- rownames(data)[distord]
  for(i in nIDs) {
    workspace.name <- paste("dist.net", i, ".RData", sep = "")
    path <- file.path(ssn.object@path, "distance", "obs",
                      workspace.name)
    if(!file.exists(path)) {
      stop("Unable to locate required distance matrix")
    }
    file_handle <- file(path, open="rb")
    distmat <- unserialize(file_handle)
    ordpi <- order(as.numeric(rownames(distmat)))
    close(file_handle)
    ni <- length(distmat[1,])
    dist.junc[(nsofar + 1):(nsofar + ni),(nsofar + 1):
                (nsofar + ni)] <- distmat[ordpi, ordpi, drop = F]
    net.zero[(nsofar + 1):(nsofar + ni),(nsofar + 1):(nsofar + ni)] <- 1
    nsofar <- nsofar + ni
  }
  
  mf <- match.call(expand.dots = FALSE)
  dataXY.out <- dataXY(formula, data,
                       family = family, trialscol = trialscol,
                       trans.power = trans.power,
                       trans.shift = trans.shift,
                       CorModels = CorModels,
                       distord = distord)
  
  REs <- dataXY.out$REs
  REmodelmatrices <- dataXY.out$REmodelmatrices
  n.all <- dataXY.out$sampsizes$n.all
  a.mat <- NULL
  b.mat <- NULL
  a.mat.data <- NULL
  b.mat.data <- NULL
  dist.hydro <- NULL
  dist.hydro.data <- NULL
  w.matrix.data <- NULL
  ind <- dataXY.out$indvecs$ind.allxy
  
  ## create any necessary matrices from distance and flow matrices
  
  if(!is.null(dist.junc) ) {
    ## maximum distance to common junction between two sites
    a.mat <- pmax(dist.junc,t(dist.junc))
    a.mat.data <- a.mat[ind,ind]
    ## minimum distance to common junction between two sites
    b.mat <- pmin(dist.junc,t(dist.junc))
    b.mat.data <- b.mat[ind,ind]
    ## hydrological distance
    dist.hydro <- as.matrix(dist.junc + t(dist.junc))
    ## subset stream distance to observed locations only
    dist.hydro.data <- dist.hydro[ind, ind]
  }
  
  if(length(grep("tailup",CorModels)) | useTailDownWeight == TRUE) {
    if(missing(addfunccol) || is.null(addfunccol) ||
       length(addfunccol) == 0 || !(addfunccol %in% colnames(data)))
      stop("The specified value for addfunccol was invalid")
    flow.con.mat <- 1 - (b.mat > 0)*1
    
    w.matrix <- sqrt(pmin(outer(data[distord,addfunccol], rep(1, times = n.all)),
                          t(outer(data[distord,addfunccol],rep(1, times = n.all)))) /
                       pmax(outer(data[distord,addfunccol],rep(1, times = n.all)),
                            t(outer(data[distord,addfunccol],rep(1, times = n.all)))))*
      flow.con.mat * net.zero
    w.matrix.data <- w.matrix[ind, ind]
  }
  net.zero.data <- net.zero[ind,ind]
  
  xcoord <- data[distord,xcol]
  ycoord <- data[distord,ycol]
  xcoord.data <- xcoord[ind]
  ycoord.data <- ycoord[ind]
  
  z <- dataXY.out$respvecs$z
  X2 <- dataXY.out$Xmats$X2
  n.allxy <- dataXY.out$sampsizes$n.allxy
  
  trialsvec <- NULL
  ## Initial parameter estimates, fixed effects, pseudo-data
  if(family == "gaussian"){
    A.5 <-  NULL
    Del.i <-  NULL
    zt <- z
  }
  
  loglik.environment <- environment()
  assign("RESULT",NULL,loglik.environment)
  
  if(length(grep("iso", tolower(CorModels))) > 0) para.length <- 7
  if(length(grep("el", tolower(CorModels))) > 0) para.length <- 6
  if(length(grep("svma", tolower(CorModels))) > 0) para.length <- 6
  
  params <- rep(0, para.length)
  TH.scale <- rep("norm", para.length)
  loglik <- 10e10
  grid_depth <- grid_depth
  theta <- rep(NA, para.length)
  
  if(para.length == 6) {
    gamma0 <- seq(lower.bound[1], upper.bound[1], length.out = grid_depth)
    range0 <- seq(lower.bound[2], upper.bound[2], length.out = grid_depth)
    phi0 <-   seq(lower.bound[3], upper.bound[3], length.out = grid_depth)
    gamma1 <- seq(lower.bound[4], upper.bound[4], length.out = grid_depth)
    range1 <- seq(lower.bound[5], upper.bound[5], length.out = grid_depth)
    phi1 <-   seq(lower.bound[6], upper.bound[6], length.out = grid_depth)
    
    for (i in 1:grid_depth) {
      for (j in 1:grid_depth) {
        for (k in 1:grid_depth) {
          for (l in 1:grid_depth) {
            for (m in 1:grid_depth) {
              for (n in 1:grid_depth) {
                theta[1] <- gamma0[i]
                theta[2] <- range0[j]
                theta[3] <- phi0[k]
                theta[4] <- gamma1[l]
                theta[5] <- range1[m]
                theta[6] <- phi1[n]
                
                loglik_temp <- try(m2LL.stream(theta=theta, m2LLdata = zt, X=X2, dist.hydro = dist.hydro.data, dichvar=dichvar,
                                               a.mat = a.mat.data, b.mat = b.mat.data, net.zero = net.zero.data,
                                               weight = w.matrix.data, x.dat = xcoord.data, y.dat = ycoord.data,
                                               Del.i = Del.i, A.5 = A.5,
                                               CorModels = CorModels, use.nugget=use.nugget,
                                               use.anisotropy=use.anisotropy, useTailDownWeight=useTailDownWeight,
                                               EstMeth = EstMeth, loglik.environment=loglik.environment,
                                               REs=NULL, scale=TH.scale, maxrang = NULL, data = data)) 
                if (loglik_temp < loglik) {
                  loglik <- loglik_temp
                  params <- c(gamma0[i], range0[j], phi0[k], gamma1[l], range1[m], phi1[n])
                }
              }
            }
          }
        }
      }
      #print(i)
    }
  }
  
  if(para.length == 7) {
    gamma0 <- seq(lower.bound[1], upper.bound[1], length.out = grid_depth)
    range0 <- seq(lower.bound[2], upper.bound[2], length.out = grid_depth)
    phi0 <-   seq(lower.bound[3], upper.bound[3], length.out = grid_depth)
    gamma1 <- seq(lower.bound[4], upper.bound[4], length.out = grid_depth)
    range1 <- seq(lower.bound[5], upper.bound[5], length.out = grid_depth)
    phi1 <-   seq(lower.bound[6], upper.bound[6], length.out = grid_depth)
    kappa <-  seq(lower.bound[7], upper.bound[7], length.out = grid_depth)
    
    for (i in 1:grid_depth) {
      for (j in 1:grid_depth) {
        for (k in 1:grid_depth) {
          for (l in 1:grid_depth) {
            for (m in 1:grid_depth) {
              for (n in 1:grid_depth) {
                for (o in 1:grid_depth) {
                  theta[1] <- gamma0[i]
                  theta[2] <- range0[j]
                  theta[3] <- phi0[k]
                  theta[4] <- gamma1[l]
                  theta[5] <- range1[m]
                  theta[6] <- phi1[n]
                  theta[7] <- kappa[o]
                  
                  loglik_temp <- try(m2LL.stream(theta=theta, m2LLdata = zt, X=X2, dist.hydro = dist.hydro.data, dichvar=dichvar,
                                                 a.mat = a.mat.data, b.mat = b.mat.data, net.zero = net.zero.data,
                                                 weight = w.matrix.data, x.dat = xcoord.data, y.dat = ycoord.data,
                                                 Del.i = Del.i, A.5 = A.5,
                                                 CorModels = CorModels, use.nugget=use.nugget,
                                                 use.anisotropy=use.anisotropy, useTailDownWeight=useTailDownWeight,
                                                 EstMeth = EstMeth, loglik.environment=loglik.environment,
                                                 REs=NULL, scale=TH.scale, maxrang = NULL, data = data)) 
                  if (loglik_temp < loglik) {
                    loglik <- loglik_temp
                    params <- c(gamma0[i], range0[j], phi0[k], gamma1[l], range1[m], phi1[n], kappa[o])
                  }
                }
              }
            }
          }
        }
      }
      #print(i)
    }
  }
    print(c(loglik, params))            
}
  


dataXY <- function(formula, data,family = "gaussian", trialscol = NULL,
                   trans.power = NULL, trans.shift = 0, distord, CorModels) {
  ##start function here, pass distord to function
  ## data.orig <- data ## not used again so commented out
  ##make sure the order of the data matches order in distance matrix
  data <- data[distord,]
  # get a list of response and covariate names
  mod.names <- as.character(attr(terms(formula, data = data),"variables"))
  tt <- model.matrix(formula,data)
  if(!all(mod.names[-1] %in% colnames(data))) stop("The specified formula causes problems in the dataXY function, try to use a simple formula, avoid indicator functions if you can.")
  ## Check random and fixed effects are different
  if(any(mod.names[-1] %in% CorModels)) stop("Random effects and fixed effects overlap")
  ## get the number of names ( + 1, as the first is always "list")
  nc.tmp <- length(mod.names)
  # name of the response variable
  response.col <- mod.names[2]
  
  # total number of observations
  n.all <- length(data[,1])
  # create a vector of all TRUE values
  ind.all <- rep(TRUE, times = n.all)
  ind.allcov <- ind.all
  # if there are any covariates ...
  if(nc.tmp > 2) {
    # create a FALSE for a record with missing values of the covariates
    for(i in 3:nc.tmp) ind.allcov <- ind.allcov & !is.na(data[,mod.names[i]])
  }
  ## Check the random effects also
  REind <- which(names(data) %in% CorModels)
  if(length(REind)) {
    for(ii in REind) ind.allcov <- ind.allcov & !is.na(data[,ii])
  }
  
  # sample size without missing covariates
  n.allcov <- sum(ind.allcov)
  # remove records that had any missing values for any of the covariates
  data1 <- data[ind.allcov,]
  #standardize numeric covariates and turn character fields to factors
  covnames <- NULL
  if(nc.tmp > 2) covnames <- mod.names[3:nc.tmp]
  #  StdXDataFrame <- NULL
  if(!is.null(covnames)) {
    for(i in 1:length(covnames)) {
      #      if(is.numeric(data1[,covnames[i]])) {
      #        xmean <- mean(data1[,covnames[i]])
      #        xsdev <- sqrt(var(data1[,covnames[i]]))
      #        data1[,covnames[i]] <- (data1[,covnames[i]] - xmean)/xsdev
      #        StdXDataFrame <- rbind(StdXDataFrame,
      #                               data.frame(VariableName = covnames[i], Mean = xmean, StdDev = xsdev))
      #      }
      if(is.character(data1[,covnames[i]])) {
        data1[,covnames[i]] <- as.factor(data1[,covnames[i]])
      }
    }
  }
  # replace response column with all 1's for design matrix of all records
  z1 <- data1[,response.col]
  data1[,response.col] <- rep(1, times = n.allcov)
  mf <- model.frame(formula, data = data1)
  mt <- attr(mf, "terms")
  X1 <- model.matrix(mt, mf)
  # get names for all terms and interactions, including those set to zero
  terms.list <- attr(mt,"term.labels")
  if(attr(mt,"intercept") == 1) effnames <- "(Intercept)"
  else effnames <- "NULL"
  if(!is.null(covnames)) {
    for (i in 1:length(terms.list)) {
      form1 <- formula(paste("~ ", terms.list[i], " - 1"))
      Xt <- model.matrix(form1, data = data1)
      effnames <- c(effnames, colnames(Xt))
    }
  }
  
  setzero <-  !(effnames %in% colnames(X1))
  # check for estimability
  X1re <- rref(X1)
  # create X1 with estimable functions
  cutX1toX2 <- apply(abs(X1re),2,sum) == 1
  if(any(cutX1toX2 == FALSE)) X1 <- X1[,cutX1toX2]
  setNA <- !(effnames %in% colnames(X1)) & !setzero
  
  data1[,response.col] <- z1
  
  #indicator vector of all records without missing covariates or response
  ind.allxy <- ind.allcov & !is.na(data[,response.col])
  #sample size of all records without missing covariates or response
  n.allxy <- sum(ind.allxy)
  # indicator vector of data1 without missing covariates
  ind.ysubx <- !is.na(data1[response.col])
  # data set of all records without missing covariates or response
  data2 <- data[ind.allxy,]
  # create a working column of the response variable
  data2[,"work.col"] <- data2[,response.col]
  # transformations on the working column
  if(!is.null(trans.power)) {
    if(family != "gaussian")
      return(print("Power transformation can only be used with gaussian family"))
    if(trans.power < 0) return(print("Power transformation must be > 0"))
    if(trans.power > .5) return(print("Power transformation must be < 0.5"))
    if(trans.power == 0) {
      if(any((data2[,response.col] + trans.shift) <= 0))
        return(print("Data must be > 0 to use log transformation"))
      data2[,"work.col"] <- log(data2[,response.col] + trans.shift)
    } else
      data2[,"work.col"] <- (data[,response.col] + trans.shift)^trans.power
  }
  
  X1 <- as.matrix(X1)
  X2 <- X1[ind.ysubx,]
  X2 <- as.matrix(X2)
  # check for estimability
  X2re <- rref(X2)
  # create X1 with estimable functions
  cutX1toX2 <- apply(abs(X2re),2,sum) == 1
  if(any(cutX1toX2 == FALSE)) X2 <- X2[,cutX1toX2]
  setNA2 <- !(effnames %in% colnames(X2))
  # data set of observed data only
  if(is.factor(data2[,"work.col"]))
    data2[,"work.col"] <- as.numeric(as.character(data2[,"work.col"]))
  if(is.character(data2[,"work.col"]))
    data2[,"work.col"] <- as.numeric(data2[,"work.col"])
  # vector of observed response variable
  z <- as.matrix(data2[,"work.col"], ncol = 1)
  attr(z,"pid") <- data2[,"pid"]
  trialsvec <- NULL
  # if response variable is binomial with n trials change z to proportion
  if(!is.null(trialscol) & family == "binomial"){
    if(is.factor(data2[,"trialscol"]))
      data2[,"trialscol"] <- as.numeric(as.character(data2[,"trialscol"]))
    if(is.character(data2[,"trialscol"]))
      data2[,"trialscol"] <- as.numeric(data2[,"trialscol"])
    trialsvec <- data2[,"trialscol"]
    z <- as.matrix(z/trialsvec, ncol = 1)
  }
  # else if response is Bernoulli, set trialsvec to all ones
  if(is.null(trialscol) & family == "binomial")
    trialsvec <- rep(1, times = n.allxy)
  
  # if any missing response values, design matrix of missing records
  data.na <- NULL
  ind.RespNA <- NULL
  if((n.allcov - n.allxy) > 0) {
    data.na <- data1[!ind.ysubx,]
    ind.RespNA <- data[,"pid"] %in% data.na[,"pid"]
  }
  # get the rank of the design matrix
  p <- sum(svd(X1)$d>1e-10)
  ## GET REs
  REs <- NULL
  REmodelmatrices <- NULL
  REind <- which(names(data) %in% CorModels)
  if(length(REind)) {
    REs <- list()
    REnames <- sort(names(data)[REind])
    ## model matrix for a RE factor
    for(ii in 1:length(REind)) REmodelmatrices[[REnames[ii]]] <- model.matrix(~data[ind.allxy,REnames[ii]] - 1)
    ## corresponding block matrix
    for(ii in 1:length(REind)) REs[[ii]] <- tcrossprod(REmodelmatrices[[ii]])
    names(REs) <- REnames
    names(REmodelmatrices) <- REnames
  }
  
  outpt <- list(
    datasets = list(
      data0 = data,
      data1 = data1,
      data2 = data2,
      data.na = data.na
    ),
    indvecs = list(
      ind.allcov = ind.allcov,
      ind.allxy = ind.allxy,
      ind.ysubx = ind.ysubx,
      ind.RespNA = ind.RespNA,
      setzero = setzero,
      setNA = setNA,
      setNA2 = setNA2,
      cutX1toX2 = cutX1toX2
    ),
    sampsizes = list(
      n.all = n.all,
      n.allcov = n.allcov,
      n.allxy = n.allxy
    ),
    Xmats = list(
      X1 = X1,
      X2 = X2,
      p = p,
      ##	  StdXDataFrame = StdXDataFrame,
      effnames = effnames
    ),
    respvecs = list(
      response.col = response.col,
      z = as.matrix(z, ncol = 1),
      trialsvec = trialsvec
    ),
    REs = REs,
    REmodelmatrices = REmodelmatrices
  )
  
  outpt
}



rref <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, fractions=FALSE) {
  
  ## A: coefficient matrix
  ## tol: tolerance for checking for 0 pivot
  ## verbose: if TRUE, print intermediate steps
  ## fractions: try to express nonintegers as rational numbers
  ## Written by John Fox
  
  if (fractions) {
    ## Added MASS to Depends in DESCRIPTION file
    ##mass <- require(MASS)
    ##if (!mass) stop("fractions=TRUE needs MASS package")
  }
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  
  n <- nrow(A)
  m <- ncol(A)
  for (i in 1:min(c(m, n))){
    col <- A[,i]
    col[1:n < i] <- 0
    # find maximum pivot in current column at or below current row
    which <- which.max(abs(col))
    pivot <- A[which, i]
    
    if (abs(pivot) <= tol) next     # check for 0 pivot
    if (which > i) A[c(i, which),] <- A[c(which, i),]  # exchange rows
    
    A[i,] <- A[i,]/pivot            # pivot
    row <- A[i,]
    A <- A - outer(A[,i], row)      # sweep
    A[i,] <- row                    # restore current row
    if (verbose)
      if (fractions) print(MASS::fractions(A))
    else print(round(A,round(abs(log(tol,10)))))
  }
  
  for (i in 1:n)
    if (max(abs(A[i,1:m])) <= tol)
      A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom
  if (fractions) fractions (A)
  else round(A, round(abs(log(tol,10))))
}



theta.ini <- function(z, X, dichvar, CorModels, use.nugget, use.anisotropy, dist.hydro.data, x.dat, y.dat, REs, init){
  n.models <- length(CorModels)
  var.resid <- mean((z - X %*% mginv(t(X) %*% X) %*% t(X) %*% z)^2)
  theta <- NULL
  scale <- NULL ## What scale is the parameter on - to simplify back transformation
  type <- NULL
  terms <- NULL
  
  if(length(grep("tailup",CorModels)) > 0){
    if(length(grep("tailup",CorModels)) > 1)
      stop("Cannot have more than 1 tailup model")
    theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                                  log(mean(dist.hydro.data/updist_scaler))),ncol = 1))
    scale <- c(scale,c("norm","log"))
    type <- c(type,c("gamma0","range0"))
    #terms <- c(terms,rep(CorModels[grep("tailup",CorModels)], times = 2))
    terms <- c(terms,c('psill', 'range'))
  }
  
  if(length(grep("taildown",CorModels)) > 0){
    if(length(grep("taildown",CorModels)) > 1)
      stop("Cannot have more than 1 taildown model")
    theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                                  log(mean(dist.hydro.data/updist_scaler))),ncol = 1))
    scale <- c(scale,c("norm","log"))
    type <- c(type,c("gamma0","range0"))
    #terms <- c(terms,rep(CorModels[grep("taildown",CorModels)], times = 2))
    terms <- c(terms,c('psill', 'range'))
  }
  
  if(length(grep("isotropic",CorModels)) > 0){
    if(length(grep("isotropic",CorModels)) > 1)
      stop("Cannot have more than 1 isotropic model")
    theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                                  log(mean(dist.hydro.data/updist_scaler))),ncol = 1))
    scale <- c(scale,c("norm","log"))
    type <- c(type,c("gamma0","range0"))
    #terms <- c(terms,rep(CorModels[grep("taildown",CorModels)], times = 2))
    terms <- c(terms,c('psill', 'range'))
  }
  
  if(use.nugget == TRUE) {
    if(is.null(theta)) {
      theta <- log(var.resid)
      scale <- "norm"
      type <- "phi0"
      terms <- "nugget"
    } else {
      theta <- c(theta,log(.1*var.resid))
      scale <- c(scale,"norm")
      type <- c(type,"phi0")
      terms <- c(terms,"nugget")
    }
  }

  #### add initial values for nonstationary params ####
  if (length(init) == 3) {
    theta <- c(theta, init)
    scale <- c(scale, rep("norm", times = 3))
    type <- c(type, c("gamma1", "range1", "phi1"))
    terms <- c(terms, "psill", "range", "nugget")
  }

  if (length(init) == 4) {
    theta <- c(theta, init)
    scale <- c(scale, rep("norm", times = 4))
    type <- c(type, c("gamma1", "range1", "phi1", "kappa"))
    terms <- c(terms, "psill", "range", "nugget", "kappa")
  }

  attr(theta,"scale") <- scale
  attr(theta,"type") <- type
  attr(theta,"terms") <- terms
  
  theta
}



mginv <- function(X, tol = sqrt(.Machine$double.eps)) { 
  dnx <- dimnames(X) 
  if(is.null(dnx)) dnx <- vector("list", 2) 
  s <- svd(X) 
  nz <- s$d > tol * s$d[1] 
  structure( 
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X, 
    dimnames = dnx[2:1]) 
}



m2LL.stream <- function(theta, m2LLdata, X, dist.hydro = NULL, dichvar=dichvar,
                        a.mat = NULL, b.mat = NULL, net.zero = NULL,
                        weight = NULL, x.dat = NULL, y.dat = NULL,
                        Del.i = NULL, A.5 = NULL,
                        CorModels, use.nugget, use.anisotropy, useTailDownWeight,
                        EstMeth = "REML",loglik.environment=loglik.environment,
                        REs=NULL, scale, maxrang = NULL, data = data) {
  if((max(theta) > 20) | (min(theta) < -20)) return(1e+32)
  theta1 <- untrans.theta(theta = theta, scale = scale)
  if(!is.null(maxrang)) 
    if (any(theta1[!is.na(maxrang)] > maxrang[!is.na(maxrang)])) return(1e+32)
  
  z <- m2LLdata
  n <- length(z)
  p <- length(X[1,])
  
  if(length(CorModels) == 0) V <- diag(rep(theta1, times = n)) 
  else V <- makeCovMat_nonstationary(theta1, dichvar = dichvar,
                       dist.hydro = dist.hydro, w.matrix = weight,
                       a.mat = a.mat, b.mat = b.mat, net.zero = net.zero,
                       x.row = x.dat, y.row = y.dat,
                       x.col = x.dat, y.col = y.dat,
                       CorModels = CorModels, useTailDownWeight = useTailDownWeight,
                       use.nugget = use.nugget, use.anisotropy = use.anisotropy, REs = REs, data = data)
  
  print(theta1)
  V.eigenvalues <- eigen(V, symmetric=TRUE, only.values=TRUE)
  
  if(any(V.eigenvalues$values <=0)) stop("covariance matrix is not positive definite")
  if(!is.null(Del.i)) V <- Del.i*A.5*t((Del.i*A.5) * V)
  
  qrV <- qr(V)
  ViX <- solve(qrV,X)
  covbi <- crossprod(X,ViX) ## Computationally more efficient than covbi <- t(X) %*% ViX
  covb <- solve(covbi)
  b.hat <- covb %*% crossprod(ViX,z) ##b.hat <- covb %*% t(ViX) %*% z
  r <- z - X %*% b.hat
  f1 <- t(r) %*% solve(qrV,r) + sum(log(abs(diag(qr.R(qrV)))))
  
  if(EstMeth == "REML") f1 <- f1 + sum(log(svd(covbi)$d))
  nmult <- (n - p)*log(2*pi)
  
  if(EstMeth == "ML") nmult <- n*log(2*pi)
  result <- f1 + nmult
  
  ## Add this value of theta and the -loglik to an environment for examination later
  
  RESULT <- get("RESULT",loglik.environment)
  RESULT <- rbind(RESULT,c(theta1,result))
  assign("RESULT",RESULT,loglik.environment)
  
  result
}



untrans.theta <- function(theta, scale) {
  ind <- which(scale == "log")
  if (length(ind)) 
    theta[ind] <- exp(theta[ind])
  ind <- which(scale == "logistic")
  if (length(ind)) 
    theta[ind] <- exp(theta[ind])/(1 + exp(theta[ind]))
  ind <- which(scale == "logistic180")
  if (length(ind)) 
    theta[ind] <- 180 * exp(theta[ind])/(1 + exp(theta[ind]))
  theta
}



makeCovMat_nonstationary <- function(theta, dichvar, dist.hydro, a.mat, b.mat, w.matrix = NULL,
                                     net.zero, x.row, y.row, x.col, y.col, useTailDownWeight,
                                     CorModels, use.nugget, use.anisotropy, REs, data = data) 
{
  nRow <- length(x.row)
  nCol <- length(x.col)
  if(is.null(net.zero)) net.zero <- matrix(1, nrow = nRow, ncol = nCol)
  V <- matrix(0, nrow = nRow, ncol = nCol)
  
  # create covariance matrix component for tailup models
  npar.sofar <- 0
  if(length(grep("tailup",CorModels)) > 0){
    if(length(grep("tailup",CorModels)) > 1)
      stop("Cannot have more than 1 tailup model")
    nonstat_type <- unlist(strsplit(CorModels,".", 
                                           fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                                 fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                                 "tailup"] - 2]
    cov_type <- substr(unlist(strsplit(CorModels,".", 
                                       fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                             fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                             "tailup"] - 1], 1, 3)
    funname <- tolower(paste(nonstat_type, ".", cov_type, ".tailup", sep = ""))
    
    tailupmod <- call(funname, dist.hydro=dist.hydro, weight = w.matrix, 
                      parsil0 = theta[npar.sofar + 1], range0 = theta[npar.sofar + 2], 
                      parsil1 = theta[npar.sofar + 4], range1 = theta[npar.sofar + 5], data = data)
    V <- V + eval(tailupmod)*net.zero
    npar.sofar <- npar.sofar + 2
  }
  # create covariance matrix component for taildown models
  if(length(grep("taildown",CorModels)) > 0){
    if(length(grep("taildown",CorModels)) > 1)
      stop("Cannot have more than 1 taildown model")
    nonstat_type <- unlist(strsplit(CorModels,".", 
                                    fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                          fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                          "taildown"] - 2]
    cov_type <- substr(unlist(strsplit(CorModels,".", 
                                       fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                             fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                             "taildown"] - 1], 1, 3)
    funname <- tolower(paste(nonstat_type, ".", cov_type, ".taildown", sep = ""))

    taildnmod <- call(funname, dist.hydro=dist.hydro, 
                      a.mat = a.mat, b.mat = b.mat, 
                      useTailDownWeight = useTailDownWeight, weight = w.matrix,
                      parsil0 = theta[npar.sofar + 1], range0 = theta[npar.sofar + 2], 
                      parsil1 = theta[npar.sofar + 4], range1 = theta[npar.sofar + 5], data = data)
    V <- V + eval(taildnmod)*net.zero
    npar.sofar <- npar.sofar + 2
  }
  # create covariance matrix component for isotropic models
  if(length(grep("isotropic",CorModels)) > 0){
    if(length(grep("isotropic",CorModels)) > 1)
      stop("Cannot have more than 1 isotropic model")
    cov_type <- substr(unlist(strsplit(CorModels,".", 
                                       fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
                                                                             fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
                                                                                             "isotropic"] - 1], 1, 3)
    funname <- tolower(paste(cov_type, ".isotropic", sep = ""))

    isomod <- call(funname, dist.hydro=dist.hydro, weight = w.matrix, 
                   parsil0 = theta[npar.sofar + 1], range0 = theta[npar.sofar + 2], 
                   parsil1 = theta[npar.sofar + 4], range1 = theta[npar.sofar + 5], kappa = theta[npar.sofar + 7], data = data)
    
    V <- V + eval(isomod)*net.zero
    npar.sofar <- npar.sofar + 2
  }
  
  # create diagonal covariance matrix component for nugget effect
  if(use.nugget == TRUE) {
    if(nRow != nCol) stop(		
      "covariance matrix asymmetric -- cannot use nugget")
    npar.sofar <- npar.sofar + 1
    V <- V + diag(exp(theta[npar.sofar] + theta[npar.sofar+3]*(data[, dichvar]/updist_scaler)), nrow = nRow, ncol = nCol)
  } else if(nRow == nCol){
    V + diag(1e-6, nrow = nRow, ncol = nCol)
  }
  
  V
}

##### Elastic model related functions ######
el.exp.tailup <- function(dist.hydro, weight, parsil0 = parsil0, range0 = range0, parsil1 = parsil1, range1 = range1, data = data)
{
  Parsil = sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Parsil*exp(-3*power.transf(ssn.assist, range1) / range0)*weight
}


el.exp.taildown <- function(dist.hydro, a.mat, b.mat, useTailDownWeight, weight = NULL, parsil0 = parsil0,
                         range0 = range0, parsil1 = parsil1, range1 = range1, data = data){
  Parsil = sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  V <- Parsil*exp(-3*power.transf(ssn.assist, range1) / range0)
  if(useTailDownWeight == TRUE) V <- V*weight
  V
}


el.lin.tailup <- function(dist.hydro, weight, parsil0 = parsil0, range0 = range0, parsil1 = parsil1, range1 = range1, data = data) {
  no <- length(dist.hydro[,1])
  np <- length(dist.hydro[1,])
  Parsil = sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Parsil*((matrix(rep(1, times = no*np), nrow = no) - power.transf(ssn.assist, range1)/range0)*(power.transf(ssn.assist, range1) < range0))*weight
}


el.lin.taildown <- function(dist.hydro, a.mat, b.mat, parsil0 = parsil0, range0 = range0, useTailDownWeight, weight = NULL,
                            parsil1 = parsil1, range1 = range1, data = data)
{
  dist.junc_law <- power.transf.tri(ssn.assist, range1)
  a_mat <- pmax(dist.junc_law, t(dist.junc_law))
  b_mat <- pmin(dist.junc_law, t(dist.junc_law))
  
  flow.connect <- b_mat == 0
  no <- length(a_mat[,1])
  np <- length(a_mat[1,])
  Parsil = sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler) )))
  
  V <- Parsil*(matrix(rep(1, times = no*np), nrow = no) - power.transf(ssn.assist, range1)/range0)*
    (power.transf(ssn.assist, range1) < range0)*flow.connect +
    Parsil*(matrix(rep(1, times = no*np), nrow = no) - a_mat/range0)*
    (a_mat < range0)*(1 - flow.connect)
  if(useTailDownWeight == TRUE) V <- V*weight
  V
}


el.sph.tailup <- function(dist.hydro, weight, parsil0 = parsil0, range0 = range0, parsil1 = parsil1, range1 = range1, data = data) {
  dist.hydro <- power.transf(ssn.assist, range1)
  no <- length(dist.hydro[,1])
  np <- length(dist.hydro[1,])
  Parsil = sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Parsil*(matrix(rep(1, times = no*np), nrow = no) - 1.5*dist.hydro/range0 + 0.5*(dist.hydro/range0)^3)*
    (dist.hydro < range0)*weight
}


el.sph.taildown <- function(dist.hydro, a.mat, b.mat, parsil0 = parsil0, range0 = range0, useTailDownWeight, weight = NULL,
                         parsil1 = parsil1, range1 = range1, data = data)
{
  dist.hydro <- power.transf(ssn.assist, range1)
  dist.junc_law <- power.transf.tri(ssn.assist, range1)
  a_mat <- pmax(dist.junc_law, t(dist.junc_law))
  b_mat <- pmin(dist.junc_law, t(dist.junc_law))
  
  flow.connect <- b_mat == 0
  no <- length(a_mat[,1])
  np <- length(a_mat[1,])
  Parsil = sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  
  V <- Parsil*(matrix(rep(1, times = no*np), nrow = no) - 1.5*dist.hydro/range0 + 0.5*(dist.hydro/range0)^3)*
    (dist.hydro < range0)*flow.connect +
    Parsil*(matrix(rep(1, times = no*np), nrow = no) - 1.5*b_mat/range0 + 0.5*a_mat/range0)*
    (matrix(rep(1, times = no*np), nrow = no) - a_mat/range0)^2*
    (a_mat < range0)*(1 - flow.connect)
  if(useTailDownWeight == TRUE) V <- V*weight
  V
}






power.transf <- function(ssn.assist, psi) {
  ssn.assistt <- ssn.assist
  if (psi != 0) {
    ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream <-
      (ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream/updist_scaler)^psi
    ssn.assistt@network.line.coords$DistanceUpstream <-
      (ssn.assistt@network.line.coords$DistanceUpstream/updist_scaler)^psi
  } else {
    ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream <-
      log(ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream/updist_scaler)
    ssn.assistt@network.line.coords$DistanceUpstream <-
      log(ssn.assistt@network.line.coords$DistanceUpstream/updist_scaler)
  }
  stream.dist.tri <- createDistMat1(ssn.assistt)
  stream.dist <- stream.dist.tri + t(stream.dist.tri)
  stream.dist <- stream.dist[-c(59,64,87,88),-c(59,64,87,88)]
  stream.dist
}



power.transf.tri <- function(ssn.assist, psi) {
  ssn.assistt <- ssn.assist
  if (psi != 0) {
    ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream <-
      (ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream/updist_scaler)^psi
    ssn.assistt@network.line.coords$DistanceUpstream <-
      (ssn.assistt@network.line.coords$DistanceUpstream/updist_scaler)^psi
  } else {
    ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream <-
      log(ssn.assistt@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream/updist_scaler)
    ssn.assistt@network.line.coords$DistanceUpstream <-
      log(ssn.assistt@network.line.coords$DistanceUpstream/updist_scaler)
  }
  stream.dist.tri <- createDistMat1(ssn.assistt)
  stream.dist.tri <- stream.dist.tri[-c(59,64,87,88),-c(59,64,87,88)]
  stream.dist.tri
}



createDistMat1 <- function (ssn) {
  driver <- RSQLite::SQLite()
  connect.name <- file.path(ssn@path, "binaryID.db")
  connect <- dbConnect(SQLite(), connect.name)
  on.exit({
    dbDisconnect(connect)
  })
  if (file.exists(file.path(ssn@path, "binaryID.db")) == FALSE)
    stop("binaryID.db is missing from ssn object")
  
  ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID <- as.factor(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID)
  net.count <- length(levels(ssn@network.line.coords$NetworkID))
  warned.overwrite <- FALSE
  for (i in 1:net.count) {
    net.num <- levels(ssn@network.line.coords$NetworkID)[i]
    ind.obs <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == as.numeric(net.num)
    site.no <- nrow(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs, ])
    
    if (site.no > 0) {
      obs.pids <- sort(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs, ])))
      net.name <- paste("net", net.num, sep = "")
      bin.table <- dbReadTable(connect, net.name)
      
      current_distance_matrix <- matrix(NA, nrow = site.no,
                                        ncol = site.no, dimnames = list(obs.pids, obs.pids))
      diag(current_distance_matrix) <- 0
      rownames(current_distance_matrix) <- obs.pids
      colnames(current_distance_matrix) <- obs.pids
      locID.obi <- attributes(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs, ])$locID
      ob.i <- as.data.frame(cbind(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,
      ])), as.numeric(levels(ssn@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[ind.obs]))[ssn@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[ind.obs]],
      locID.obi[ind.obs]))
      colnames(ob.i) <- c("pid", "rid", "locID")
      ob.i$locID <- as.factor(ob.i$locID)
      ob.i$binaryID <- bin.table$binaryID[match(ob.i$rid,
                                                bin.table$rid)]
      ob.i <- ob.i[order(ob.i[, "pid"]), ]
      rownames(ob.i) <- ob.i$pid
      ob.i_by_locID <- ob.i[order(ob.i[, "locID"]), ]
      ob.i_by_locID$pid <- as.numeric(ob.i_by_locID$pid)
      ob.i_by_locID$locID <- as.numeric(ob.i_by_locID$locID)
      ob.j_reordering <- order(ob.i_by_locID$pid)
      locID.old <- -1
      ind.dup <- !duplicated(ob.i_by_locID$locID)
      for (j in 1:nrow(ob.i)) {
        pid.i <- ob.i[j, "pid"]
        locID.i <- ob.i[j, "locID"]
        if (locID.i != locID.old) {
          junk <- get.rid.fc(ob.i_by_locID[ind.dup, "binaryID"],
                             ob.i$binaryID[j])
          ob.j <- getObsRelationshipsDF(ssn, pid.i, junk,
                                        ind.dup, ob.i, ob.i_by_locID, bin.table)
          upDist.i <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[paste(pid.i), "DistanceUpstream"]
          ob.j <- ob.j[ob.j_reordering, ]
          ind.fc <- ob.j$fc == 1
          dist.obs <- ifelse(ind.fc, upDist.i - ob.j$upDist.j,
                             upDist.i - ob.j$juncDist)
          current_distance_matrix[, paste(pid.i)] <- ifelse(dist.obs < 0, 0, dist.obs)
        }
        else {
          current_distance_matrix[, paste(pid.i)] <- current_distance_matrix[,  paste(pid.old)]
        }
        pid.old <- pid.i
        locID.old <- locID.i
      }
    }
  }
  current_distance_matrix
}



getObsRelationshipsDF <- function(ssn, pid, junk, ind, ob, ob_by_locID, bin) {
  
  ## ssn = SpatialStreamNetwork object
  ## pid = integer pid value of interest
  ## junk = ox2 data.frame with columns fc (logical) and binaryID
  ## ind = vector indicator saying whether it is a duplicate
  ## ob = data.frame with pid, rid, locID, and binaryID for sites on network
  ##      ordered by pid
  ## ob_by_locID = data.frame with pid, rid, locID, and binaryID for sites on network
  ##               ordered by locID
  ## bin = binaryID table
  
  ## Returns a data.frame that relates all sites to pid.i:
  ## pid: numeric
  ## locID: numeric
  ## fc: logical - is the sute fc with the pid of interest
  ## binaryID: binaryID of the common downstream junction
  ## junc.rid: rid for the common downstream junction
  ## upDist.j: upDist for each site
  
  ## Create relationships table
  
  ob.j.r <- data.frame(ob_by_locID[ind, c("pid", "locID")], junk,
                       stringsAsFactors = FALSE)
  
  ob.j.r$fc <- as.logical(ob.j.r$fc)
  rownames(ob.j.r)<- ob.j.r$pid
  
  ## Add column showing rid for common downstream junction (ob.j.r relates all               ##sites to pid)
  ob.j.r$junc.rid <- bin$rid[match(ob.j.r$binaryID, bin$binaryID)]
  
  ## This doesn't make sense to me...
  reps <- as.numeric(ob_by_locID$locID)
  ob.j <- ob.j.r[reps,]
  
  ## Create some funky rownames, with extension .fc
  rownames(ob.j) <- paste(rownames(ob), ".fc", sep = "")
  
  ## Don't know why we're doing this...
  ob.j$pid <- ob_by_locID$pid
  
  ## juncDist is the upstream distance of the common downstream rid junction
  ob.j$juncDist <- ssn@network.line.coords$DistanceUpstream[match(ob.j$junc.rid, ssn@network.line.coords$SegmentID)]
  
  ## upDist.j is the upDist for each observed site
  ob.j$upDist.j <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream[
    match(ob.j$pid, as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords)))]
  
  ob.j
}



get.rid.fc <- function(binIDs, referenceBinID) {
  ind.match <- .Call("test_fc", binIDs, referenceBinID)
  data.frame(fc = ind.match<0,
             binaryID = substr(binIDs, 1, abs(ind.match)),
             stringsAsFactors = FALSE)
}


##### SVMA model related functions ######
svma.exp.tailup <- function(dist.hydro, weight, parsil0 = parsil0, range0 = range0, 
                       parsil1 = parsil1, range1 = range1, data = data){
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Range <- (range0 + range1 * min_updist(data[, dichvar])/updist_scaler)
  Parsil*exp(-3 * dist.hydro/updist_scaler / Range)*weight
}


svma.exp.taildown <- function(dist.hydro, a.mat, b.mat, useTailDownWeight, weight = NULL, parsil0 = parsil0, 
                         range0 = range0, parsil1 = parsil1, range1 = range1, data = data){
  flow.connect <- b.mat == 0
  
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Range_1 <- (range0 + range1 * min_updist(data[, dichvar])/updist_scaler)
  Range_2 <- (range0 + range1 * max_updist(data[, dichvar])/updist_scaler)
  
  V <- Parsil*exp(-3*dist.hydro/updist_scaler / Range_2)*flow.connect + 
    Parsil*exp(-3*(b.mat*Range_2+a.mat*Range_1)/updist_scaler / (Range_1*Range_2))*(1 - flow.connect)
}


svma.lin.tailup <- function(dist.hydro, weight, parsil0 = parsil0, range0 = range0, 
                       parsil1 = parsil1, range1 = range1, data = data) {
  no <- length(dist.hydro[,1])
  np <- length(dist.hydro[1,])
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Range_1 <- (range0 + range1 * min_updist(data[, dichvar])/updist_scaler)
  Range_2 <- (range0 + range1 * max_updist(data[, dichvar])/updist_scaler)
  
  V <- Parsil*((matrix(rep(1, times = no*np), nrow = no) - dist.hydro/updist_scaler/Range_1)*
                  (dist.hydro/updist_scaler < Range_1))
  V*weight
}


svma.lin.taildown <- function(dist.hydro, a.mat, b.mat, parsil0 = parsil0, range0 = range0, useTailDownWeight, weight = NULL,
                         parsil1 = parsil1, range1 = range1, data = data) {
  flow.connect <- b.mat == 0
  no <- length(a.mat[,1])
  np <- length(a.mat[1,])
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Range_1 <- (range0 + range1 * min_updist(data[, dichvar])/updist_scaler)
  Range_2 <- (range0 + range1 * max_updist(data[, dichvar])/updist_scaler)
  
  V <- Parsil*(matrix(rep(1, times = no*np), nrow = no) - dist.hydro/updist_scaler / Range_2)*
    (dist.hydro/updist_scaler < Range_2)*flow.connect + 
    Parsil*(matrix(rep(1, times = no*np), nrow = no) - a.mat/updist_scaler / Range_2)*
    (a.mat/updist_scaler < Range_2)*(1 - flow.connect)
  V
}


svma.sph.tailup <- function(dist.hydro, weight, parsil0 = parsil0, range0 = range0, parsil1 = parsil1, range1 = range1, data = data) {
  no <- length(dist.hydro[,1])
  np <- length(dist.hydro[1,])
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Range_1 <- (range0 + range1 * min_updist(data[, dichvar])/updist_scaler)
  Range_2 <- (range0 + range1 * max_updist(data[, dichvar])/updist_scaler)
  
  V <- Parsil*(dist.hydro/updist_scaler-Range_1+3*Range_2)/(2*Range_2)*
    ((matrix(rep(1, times = no*np), nrow = no) - dist.hydro/updist_scaler/Range_1)^2*
       (dist.hydro/updist_scaler < Range_1))
  V*weight
}


svma.sph.taildown <- function(dist.hydro, a.mat, b.mat, parsil0 = parsil0, range0 = range0, useTailDownWeight, weight = NULL,
                         parsil1 = parsil1, range1 = range1, data = data) {
  flow.connect <- b.mat == 0
  no <- length(a.mat[,1])
  np <- length(a.mat[1,])
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Range_1 <- (range0 + range1 * min_updist(data[, dichvar])/updist_scaler)
  Range_2 <- (range0 + range1 * max_updist(data[, dichvar])/updist_scaler)
  
  V <- Parsil*(dist.hydro/updist_scaler-Range_2+3*Range_1)/(2*Range_1)*
    (matrix(rep(1, times = no*np), nrow = no) - dist.hydro/updist_scaler / Range_2)^2*
    (dist.hydro/updist_scaler < Range_2)*flow.connect + 
    Parsil*(a.mat/updist_scaler-3*b.mat/updist_scaler-Range_2+3*Range_1)/(2*Range_1)*
    (matrix(rep(1, times = no*np), nrow = no) - a.mat/updist_scaler / Range_2)^2*
    (a.mat/updist_scaler < Range_2)*(1 - flow.connect)
  V
}


min_updist <- function (updist) {
  l <- length(updist)
  s1_mat <- matrix(0, l, l)
  for (i in 1:l) {
    for (j in 1:l) {
      s1_mat[i,j] <- ifelse(updist[i] < updist[j], updist[i], updist[j])
    }
  }
  s1_mat
}


max_updist <- function (updist) {
  l <- length(updist)
  s2_mat <- matrix(0, l, l)
  for (i in 1:l) {
    for (j in 1:l) {
      s2_mat[i,j] <- ifelse(updist[i] >= updist[j], updist[i], updist[j])
    }
  }
  s2_mat
}


##### Isotropic model related functions ######
exp.isotropic <- function(dist.hydro, a.mat, b.mat, useTailDownWeight, weight = NULL, parsil0 = parsil0,
                    range0 = range0, parsil1 = parsil1, range1 = range1, kappa = kappa, data = data){
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Parsil*exp(-3*(power.transf(ssn.assist, range1))^kappa / range0)
}

mat.isotropic <- function(dist.hydro, a.mat, b.mat, useTailDownWeight, weight = NULL, parsil0 = parsil0,
                    range0 = range0, parsil1 = parsil1, range1 = range1, kappa = kappa, data = data){
  f <- function(x) matern(x, 1, kappa) - 0.05
  eff_range0 <- uniroot(f, c(0,5))$root
  
  noweight <- weight
  Parsil <- sqrt(tcrossprod(exp(parsil0 + parsil1*(data[, dichvar]/updist_scaler))))
  Parsil*matern(power.transf(ssn.assist, range1), range0/eff_range0, kappa)
}