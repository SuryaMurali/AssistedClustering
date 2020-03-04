MDistance <- function(n,D)
{  #------------------------------------------------------------#
  # function to make full distance Matrix from  distance vector
  #------------------------------------------------------------#
  MD <- matrix(0, nrow=n, ncol=n)
  for(j in 1:(n-1) )
    for(k in (j+1):n)
    {    kk <- n*(j-1) - j*(j-1)/2 + k-j
    MD[j,k] <- MD[k,j] <- D[kk]
    }
  return(MD)
}    
#

MDistanceVectorFtn <- function(DistanceMatrix)
{ #---------------------------------------------------#
  # function of making distance vector(upper) from 
  # full distance matrix
  #---------------------------------------------------#
  n <- ncol(DistanceMatrix)
  mv <- array(0, n*(n-1)/2 )
  for(i in seq(n-1))
    for(j in seq(from = i+1, to = n))
    {   arr.n <- n*(i-1) - i*(i-1)/2 + j-i
    mv[arr.n] <- DistanceMatrix[i,j]
    }
  return(mv)
}


Get.Subdata <- function(data)
{
  fullsize <- nrow(data)
  nvar <- ncol(data)
  cat("\n Cases : ",fullsize,",", "Var's : ",nvar) 
  x <- 1:fullsize
  cat("\n Type Sample percent (10-100%)  : ")
  percent <- readline()
  percent <- as.integer(percent)
  size <- floor(fullsize*percent/100)
  subdata <- sample(x, size)
  subdata <- data[subdata,]
  return(subdata)
}


DecideClusterCall.ward <- function(n, DistA,wk)
{
  
  # Ward Cluster agglomeration
  ward.res <- hclust(DistA, method="ward.D")
  
  # Mojena's stopping rule 
  # cut height = mean(heights of cluster N-1) + wk*Sh(standard deviation of heights)
  # Mojena(1980) suggested alpha to be under 2.5 - Stopping rules for Ward's clustering method, 
  # COMPSTAT 1980 Proceedings, Physica-Verlag, 426-432.
  
  
  Height <- ward.res$height
  meanHeight <-  mean(Height)
  stdevHeight <- sd(Height)   # s-plus : stdev 
  cutHeight <- meanHeight + wk*stdevHeight
  
  groupNumber <- 1
  
  for(i in 1:(n-1) )
  {
    if(Height[n-i] > cutHeight) groupNumber <- groupNumber + 1
    else break
  }            
  
  # cat("\n Group Number = ", groupNumber)
  
  clsInfo <- cutree(ward.res, groupNumber)
  
  # cat("\n --- Height ---\n")
  # print(Height)
  
  # cat("\n --- cutHeight = ", cutHeight, " Group = ", groupNumber)
  
  # cat("\n --- Cluster Using Ward's Method  \n")
  # print(clsInfo)
  
  clsResult <- list(clsInfo = clsInfo, groupNumber = groupNumber)
  clsResult
}

WardClustering <- function(data,wk)
{
  # Initialization : number of clusters, cluster centroids 
  # Using Ward's method
  
  # Make Euclidean distances between observations for MST
  DistanceArray <- dist( data, method="euclidean")
  
  data = as.matrix(data) 
  n <- nrow(data)
  
  clsResult <- DecideClusterCall.ward(n, DistanceArray,wk)
  clsInfo   <- clsResult$clsInfo
  clstGroup <- clsResult$groupNumber
  
  Gmean <- matrix(0, clstGroup, ncol(data))
  
  clsdata <- cbind(clsInfo, data)
  
  
  Glength <- array(0, clstGroup)
  
  for(i in 1:clstGroup) { Glength[i] <- length(clsInfo[clsInfo == i]) }
  
  for(i in 1:clstGroup)
  {    gdata <- clsdata[clsdata[,1] == i,]
  if( Glength[i] == 1 )
  {  Gmean[i,] <- gdata[-1] }
  else
  { 
    Gm <- apply(gdata, 2, 'mean')
    Gmean[i,] <- Gm[-1] 
  }
  }
  
  Mresult <- list(clsInfo=clsInfo, clsGroup=clstGroup, clsCenter=Gmean)
  Mresult
  
}

kmeans.center <- function(clsInfo, clst.num, kdata)
{ 
  totdata <- cbind(clsInfo, kdata)
  Glength <- array(0, clst.num)
  Gmean <- matrix(0, clst.num, ncol(kdata))
  
  for(i in 1:clst.num) { Glength[i] <- length(clsInfo[clsInfo == i]) }
  for(i in 1:clst.num)
  {    gdata <- totdata[totdata[,1] == i,]
  if( Glength[i] == 1 )
  { Gmean[i,] <- gdata[-1] }
  else
  { Gm <- apply(gdata, 2, 'mean')
  Gmean[i,] <- Gm[-1] }
  }
  
  return(Gmean)
}



sample.data <- function(data)
{
  fullsize = nrow(data)
  nvar <- ncol(data)
  x <- 1:fullsize
  
  cat("\n Type Sampling Rate (10-100%, def=10%)  : ")
  percent <- readline()
  if(percent == "") percent <- 10
  percent <- as.integer(percent)
  size <- floor(fullsize*percent/100)
  
  sub <- srswor(size, length(x))
  samplelist <- (1:length(x))[sub==1]
  subdata <- data[samplelist,]
  
  return(subdata)
}

TwoStepKmeans <- function(sam.data, k.data, wk)
{
  InitRes <- WardClustering(sam.data,wk)
  clsGroup <- InitRes$clsGroup
  initCenter <- InitRes$clsCenter
  
  kmeans.out <- kmeans(k.data, initCenter)
  
  # print(kmeans.out)
  return(kmeans.out)
}

Between.distance <- function(x, center)
{
  betw.distance = sqrt(rowSums((x-center)^2 ))
  maxbdist = max(betw.distance)
  decide.ratio = 0.95
  dist.ratio = betw.distance / maxbdist
  large.distance = betw.distance[dist.ratio >= decide.ratio]
  cutoff = min(large.distance)
  bd.out = list(bd=betw.distance, cutoff = cutoff)
  return(bd.out) 
}

Mahal.distance <- function(nofx, x, quantile, cls.ratio, center)
{
  nx = nrow(x)
  dist.method = 1
  md.df = ncol(x)
  
  # delete variable which is nearly constant 
  x.cov = cov(x)
  x.diag = diag(x.cov)
  x.delete = which(x.diag <= 0.000001)
  x.length = length(x.delete)
  if(x.length > 0) x = x[,-x.delete]
  
  if(length(ncol(x)) == 0 ) ncol.x = 1
  else
    ncol.x = ncol(x)
  
  
  if( (nx / nofx) < cls.ratio ) 
  {  mdist = rep(97, times=nx)  
  dist.method = 7
  }
  else if ( nx <= ncol.x ) 
  {  mdist = rep(95, times=nx)
  dist.method = 5
  }
  else
  {
    mdist = mahalanobis(x, colMeans(x), cov(x) )
    mdist = sqrt(mdist)
    dist.method = 1
  }
  
  cutoff = qchisq(quantile, md.df)
  cutoff = sqrt(cutoff)
  
  md.out <- list(mdist=mdist, cutoff = cutoff, method=dist.method)
  return(md.out)
}

RMahal.distance <- function (X, quantile, rtype.sel) 
{
  require(robustbase)
  require(MASS)
  dmethod = 2
  md.df = ncol(X)
  
  # Robust scale parameter check !
  
  # delete variable which is nearly constant 
  X.cov = cov(X)
  X.diag = diag(X.cov)
  X.delete = which(X.diag <= 0.000001)
  X.length = length(X.delete)
  if(X.length > 0) X = X[,-X.delete]
  
  if(length(ncol(X)) == 0 )
  {
    nc = 1
    newX = X
  }
  else
  {
    nc = ncol(X)
    iqr.x = array(0, nc)
    
    for(i in 1:nc)
    {  iqr.x[i] = IQR(X[,i]) }
    
    iqr.delete = which(iqr.x == 0)
    iqr0.len = length(iqr.delete)
    
    if(iqr0.len > 0) newX = X[,-iqr.delete]
    else 
      newX = X
    
    if(length(ncol(newX)) == 0 ) nc = 1
  }
  
  
  if(nc == 1)     # check if only one-variable is remained
  {  rd = mahalanobis(X, colMeans(X), cov(X) )
  rd = sqrt(rd)  
  dmethod = 3 
  cutoff = sqrt(qchisq(quantile, md.df))
  }
  else
  { if(rtype.sel == 1) X.mcd = cov.mve(newX)
  else if (rtype.sel == 2) X.mcd = cov.rob(newX)
  else if (rtype.sel == 3) X.mcd = cov.mcd(newX)
  else
    X.mcd <- cov.Mcd(newX, 0.6)
  
  rd = sqrt(mahalanobis(newX, X.mcd$center, X.mcd$cov))
  cutoff <- sqrt(qchisq(quantile, md.df))
  }
  
  md.out = list(rd = rd, cutoff = cutoff, dmethod=dmethod)
  return(md.out)
}

RobustMahal.distance <- function(nofx, x, quantile, cls.outlier, rtype.sel, center)
{
  nx = nrow(x)
  dist.method = 1
  md.df = ncol(x)
  
  if( (nx / nofx) < cls.outlier ) 
  {   rmdist = rep(97, times=nx)  
  dist.method = 7
  }
  else if ( nx <= (ncol(x)+1) ) 
  {   rmdist = rep(95, times=nx)
  dist.method = 5
  } 
  else
  {
    rd.dist = RMahal.distance(x, quantile, rtype.sel)
    rmdist = rd.dist$rd
    dist.method=rd.dist$dmethod
  }   
  
  
  cutoff = qchisq(quantile, md.df)
  cutoff = sqrt(cutoff)
  
  md.out <- list(mdist=rmdist, cutoff = cutoff, method=dist.method)
  return(md.out)
}

Outlier.detection.fn <- function(k.data, cls.center, cls.id, cls.num, id.label, q.value, cls.ratio, mahal.select, rtype.sel)
{
  
  # Local Outlier based on Clustering Results 
  # Calculate Mahal distance for each cluster
  # cat("\n\n === Outlier Detection === \n")
  
  Ldistance = array()
  Loutlier = array()
  Ldtype = array()       # distance type 
  
  nofx = nrow(k.data)
  k1 = 1
  
  for (i in 1:cls.num)
  {   
    cls.data <- k.data[cls.id==i, ]
    rownames(cls.data) <- id.label[cls.id == i]
    rowid <- rownames(cls.data) 
    rowid = as.numeric(rowid)
    center = cls.center[i,]
    
    if(mahal.select == 1) 
      md.out = Mahal.distance(nofx, cls.data, q.value, cls.ratio, center)
    else
      md.out = RobustMahal.distance(nofx, cls.data, q.value, cls.ratio, rtype.sel, center)
    
    # md.out = Moutlier(cls.data, quantile=0.99, plot=FALSE)
    
    mdist = md.out$mdist
    Ldmethod = md.out$method
    
    # dmethod : 1=Mahal distance, 2=Robust Mahal distance, 5,7=small size cluster  
    # identify outliers 
    
    cutoff.value <- md.out$cutoff 
    
    mahal.outlier <- which(mdist >= cutoff.value)
    k = length(mahal.outlier)
    if( k > 0 )
    { k2 <- k1 + (k-1)
    Ldistance[k1:k2] =mdist[mahal.outlier]
    Loutlier[k1:k2] = rowid[mahal.outlier]
    Ldtype[k1:k2] = Ldmethod
    k1 <- k2+1 
    }  
    
    #if(Ldmethod == 9) Ldistance[Ldistance >0 ] = 99
    
  }
  
  Ldistance <- as.numeric(Ldistance)  
  Loutlier <- as.numeric(Loutlier)
  Ldistance = round(Ldistance,3)
  #Ldistance[Ldistance==99] = NA
  cutoff.value = round(cutoff.value,3)
  
  # Global Outlier based on total data sets  
  # Calculate Mahal distance for total data sets
  
  Gdistance = array()
  Goutlier = array()
  Gdtype = array()
  nofx = nrow(k.data)
  k1 = 1
  center = colMeans(k.data)
  
  if(mahal.select == 1) 
    Gmd.out = Mahal.distance(nofx, k.data, q.value, cls.ratio, center)
  else
    Gmd.out = RobustMahal.distance(nofx, k.data, q.value, cls.ratio, rtype.sel, center)
  
  Gmdist = Gmd.out$mdist
  Gdmethod = Gmd.out$method
  # identify outliers 
  
  Gcutoff.value <- Gmd.out$cutoff 
  
  Gmahal.outlier <- which(Gmdist >= Gcutoff.value)
  rowid = id.label
  
  k = length(Gmahal.outlier)
  if( k > 0 )
  { k2 <- k1 + (k-1)
  Gdistance[k1:k2] =Gmdist[Gmahal.outlier]
  Goutlier[k1:k2] = rowid[Gmahal.outlier]
  Gdtype[k1:k2] = Gdmethod
  k1 <- k2+1 
  }  
  
  #if(Gdmethod == 9) Gdistance[Gdistance >0 ] = 99
  Gdistance <- as.numeric(Gdistance)  
  Goutlier <- as.numeric(Goutlier)
  Gdistance = round(Gdistance,3)
  Gcutoff.value = round(Gcutoff.value,3)
  
  outlier.data <- list(Ldistance=Ldistance, Loutlier=Loutlier, Ldtype=Ldtype, Gdistance=Gdistance, Goutlier=Goutlier, Gdtype=Gdtype, cutoff=cutoff.value)
  return(outlier.data)
  
}

PrintOutlier.ftn <- function(Loutlier, Ldistance,Ldtype, Goutlier, Gdistance, Gdtype, cutoff, ss.ratio, mahal.select)
{
  # Local Outliers
  
  PLdtype = Ldtype
  PGdtype = Gdtype
  CLdtype = Ldtype
  CGdtype = Gdtype
  PLdtype[PLdtype==1] = "-"
  PLdtype[PLdtype==2] = "-"
  PLdtype[PLdtype==3] = "md"
  PLdtype[PLdtype==5] = "sp"
  PLdtype[PLdtype==7] = "sn"
  PGdtype[PGdtype==1] = "-"
  PGdtype[PGdtype==2] = "-"
  PGdtype[PGdtype==3] = "md"
  PGdtype[PGdtype==5] = "sp"
  PGdtype[PGdtype==7] = "sn"
  CLdtype[PLdtype==2] = 1
  CGdtype[PGdtype==2] = 1
  uniq.CL = unique(CLdtype)
  uniq.GL = unique(CGdtype)
  
  
  Lorder = order(Ldistance, decreasing=T)
  PLoutlier = Loutlier[Lorder]
  PLdistance = Ldistance[Lorder]
  PPLdtype = PLdtype[Lorder]
  Gorder = order(Gdistance, decreasing=T)
  PGoutlier = Goutlier[Gorder]
  PGdistance = Gdistance[Gorder]
  PPGdtype = PGdtype[Gorder]
  
  
  if( is.na(Loutlier[1]) == FALSE ) n.outlier = length(Loutlier) 
  else n.outlier = 0 
  
  if( n.outlier == 0)
  {
    cat("\n        Potential Outliers(Local) = - ")
    if(mahal.select == 1) cat("\n        Mahal. Distance(Local) = - ")
    else cat("\n        Robust Mahal. Distance(Local) = - ") 
    cat("\n        Number of Outliers(Local) = 0 ")
  }
  else
  {
    cat("\n        Potential Outliers(Local) = ", PLoutlier)
    if(mahal.select == 1) cat("\n        Mahal. Distance(Local) = ", PLdistance)
    else cat("\n        Robust Mahal. Distance(Local) = ", PLdistance)
    if(any(uniq.CL > 2)) cat("\n        (*)Outlier Type = ", PPLdtype)
    cat("\n        Number of Outliers(Local) = ", n.outlier)
  }
  
  
  
  
  # Global Outliers
  
  if( is.na(Goutlier[1]) == FALSE ) n.goutlier = length(Goutlier) 
  else n.goutlier = 0 
  
  if( n.goutlier == 0)
  {
    cat("\n        Potential Outliers(Global) = - ")
    if(mahal.select == 1) cat("\n        Mahal. Distance(Global) = - ")
    else cat("\n        Robust Mahal. Distance(Global) = - ") 
    cat("\n        Number of Outliers(Global) = 0 ")
  }
  else
  {
    cat("\n        Potential Outliers(Global) = ", PGoutlier)
    if(mahal.select == 1) cat("\n        Mahal. Distance(Global) = ", PGdistance)
    else cat("\n        Robust Mahal. Distance(Global) = ", PGdistance)
    if(any(uniq.GL > 2)) cat("\n        (*)Outlier Type = ", PPGdtype)
    cat("\n        Outlier Cutoff = ", cutoff )
    cat("\n        Number of Outliers(Global) = ", n.goutlier)
  }
  cat("\n        SSB/SST    = ", round(ss.ratio,4) , "\n")
  
}

RunStep5.Step6.rtn <- function(sam.data, kdata, selected.vars, unselected.vars,wk, id.label,q.value, cls.ratio, mahal.select, rtype.sel)
{
  
  no.col <- ncol(kdata)
  pselected.vars = sort(selected.vars)
  e.data = sam.data[,c(pselected.vars)]
  k.data = kdata[,c(pselected.vars)]
  kmeans.out = TwoStepKmeans(e.data, k.data,wk)
  cls.id = kmeans.out$cluster
  cls.size = kmeans.out$size
  cls.num = length(cls.size)
  cls.center = kmeans.out$centers
  
  cat("\n        Selected Var's = (", selected.vars,")" )
  unselected.num = length(unselected.vars)
  if(unselected.num > 0 ) cat("\n        UnSelected Vars = (", unselected.vars, ")" )
  else  cat("\n        UnSelected Vars = ( - )" )  
  cat("\n        Number of Cluster = ", cls.num)
  cat("\n        Cluster Sizes = ", cls.size)
  
  # Outlier Detection Procedure
  out.data = Outlier.detection.fn(k.data, cls.center, cls.id, cls.num, id.label, q.value, cls.ratio, mahal.select,rtype.sel)
  
  Loutlier = out.data$Loutlier
  Ldistance= out.data$Ldistance
  Ldtype = out.data$Ldtype
  Goutlier = out.data$Goutlier
  Gdistance= out.data$Gdistance
  Gdtype = out.data$Gdtype
  cutoff = out.data$cutoff
  ss.ratio = kmeans.out$betweenss / kmeans.out$totss
  
  PrintOutlier.ftn(Loutlier, Ldistance, Ldtype, Goutlier, Gdistance, Gdtype, cutoff, ss.ratio, mahal.select) 
  
  YPj = kmeans.out$cluster
  
  y.u.adj = array(0, no.col)
  
  for(j in c(unselected.vars) )
  {
    num.unique = length(unique(kdata[,j]))
    if( num.unique < cls.num ) Ugroup = kdata[,j]
    else
    { kmeans.v <- kmeans(kdata[,j], cls.num)
    Ugroup <- kmeans.v$cluster           
    } 
    y.u.adj[j] = adjustedRandIndex(YPj, Ugroup)  
  }
  
  
  u.adj.max = max(y.u.adj)
  u.adj.which = which.max(y.u.adj)
  
  if( u.adj.max > 0 ) cat("\n\n  Repeat(Step5-6) : adj.max = ", round(u.adj.max,4), " which =(", u.adj.which,")" )
  else cat("\n\n  Repeat(Step5-6) : adj.max = ", round(u.adj.max,4), " which = ( - ) ") 
  
  result <- list(adjmax = u.adj.max, adjwhich=u.adj.which, cls.id=cls.id, cls.size=cls.size, Loutlier=Loutlier, Ldistance=Ldistance, Ldtype=Ldtype, Goutlier=Goutlier, Gdistance=Gdistance, Gdtype=Gdtype, cutoff=cutoff, ratio=ss.ratio)
  return(result)
}

AdjustedSim.Calculate <- function(kdata, id.label, cls.size, cls.id, all.vars, selected.vars, unselected.vars, Loutlier.potential, Goutlier.potential) 
{
  #---------------------------------------------------------------------------------#
  #==== Adjusted Rand Index between input data and Kmeans result ===                #
  #---------------------------------------------------------------------------------#
  
  clusterGroup = cls.id
  cls.num = length(cls.size)
  
  # Simulated Group File Name 
  cat("\n --- original group member file : ")
  filename = readline()
  mem.name = paste(filename, ".mem", sep="")
  originGroup = scan(mem.name)
  noisy.file = paste(filename, ".noisy", sep="")
  noisy.vars = scan(noisy.file) 
  true.vars = all.vars[-noisy.vars]
  clu.size = table(originGroup)
  origin.outlier = id.label[originGroup == 0]
  sort.Loutlier = sort(Loutlier.potential, decreasing=T)
  sort.Goutlier = sort(Goutlier.potential, decreasing=T)
  
  cat("\n ----------- < Simulation Data by clusterGeneration > ------------")
  cat("\n Original True Var's = (", true.vars,")")
  cat("\n Original Noisy Var's = (", noisy.vars,")")
  cat("\n Original Cluster Size = ", clu.size)  
  cat("\n Original Outliers = ", origin.outlier) 
  
  selected.vars = sort(selected.vars)
  cat("\n\n ---------- < Automated K-Means Clustering > ------------")
  cat("\n Selected Var's = (", selected.vars,")" )
  cat("\n UnSelected Var's = (", unselected.vars, ") " )
  cat("\n Cluster Size = ", cls.size)  
  cat("\n Potential Outliers(Local) = ", sort.Loutlier) 
  cat("\n Number of Potential Outliers(Local) = ", length(Loutlier.potential) ) 
  cat("\n Potential Outliers(Global) = ", sort.Goutlier) 
  cat("\n Number of Potential Outliers(Global) = ", length(Goutlier.potential) ) 
  
  # Adjusted Rand Index between Simulated Data and Kmeans Result
  
  #----------------------------------------------------------#
  # K-means with selected variables and without outliers     #
  #----------------------------------------------------------#
  
  #Goutlier.num = length(Goutlier.potential)
  #no.row = nrow(kdata)
  
  #if(Goutlier.num > 0.3 * no.row) toutlier.potential = Loutlier.potential
  #else 
  # {  toutlier.potential = append(Loutlier.potential, Goutlier.potential)
  #    toutlier.potential = unique(toutlier.potential)
  # }
  
  toutlier.potential = Loutlier.potential
  
  if( is.na(toutlier.potential[1]) == FALSE )  
  {
    originGroup[toutlier.potential] = 0
    origin.cluster = originGroup[originGroup > 0]
    origin.id = id.label[originGroup==0]
    cluster.id = clusterGroup[originGroup > 0]
    #cat("\n origin.outlier = ", length(origin.cluster))
    #cat("\n cluster.outlier= ", length(cluster.id))
  }
  else
  {
    origin.cluster = originGroup[originGroup > 0]
    cluster.id = clusterGroup[originGroup > 0]
    #cat("\n origin.outlier2 = ", length(origin.cluster))
    #cat("\n cluster.outlier2= ", length(cluster.id))
  }
  
  adjRand2 = adjustedRandIndex(cluster.id, origin.cluster) 
  adjRand2 = round(adjRand2,4)
  
  cls.table = table(origin.cluster, cluster.id)
  
  cat("\n Adjusted Rand Index = ", adjRand2,"\n")
  cat("\n Confusion Matrix \n") 
  print(cls.table)
}

# sum of squares
ss.value <- function(x) sum(scale(x, scale = FALSE)^2)

ratio.calculate <- function(zdata, mem)
{
  clst.num = length(unique(mem)) 
  col.num = ncol(zdata)
  row.num = nrow(zdata)
  fitted.x = matrix(0, ncol=col.num, nrow=row.num)
  zdata.center = kmeans.center(mem, clst.num, zdata)
  #print(zdata.center)
  
  for(i in 1:row.num) fitted.x[i,] = zdata.center[mem[i],] 
  fitted.x = data.frame(fitted.x)
  
  resid.x <- zdata - fitted.x
  
  #print(head(fitted.x))
  
  ratio.val = ss(fitted.x) / ss(zdata)
  #cat("\n Total SS   = ", ss.value(zdata) , "\n")
  #cat("\n Between SS = ", ss.value(fitted.x) , "\n")
  #cat("\n Within  SS = ", ss.value(resid.x) , "\n")
  #cat("\n SSB/SST    = ", round(ratio.val,3) , "\n")
  return(ratio.val)
}



AdjustedRdata.Calculate <- function(kdata, cls.size, cls.id, selected.vars, Loutlier.potential,Goutlier.potentialG) 
{
  #---------------------------------------------------------------------------------#
  #==== Adjusted Rand Index between input data and Kmeans result ===                #
  #---------------------------------------------------------------------------------#
  
  clusterGroup = cls.id
  cls.num = length(cls.size)
  
  # Real Data
  cat("\n --- Original Group Member file : ")
  mem.name = readline()
  originGroup = scan(mem.name)
  
  cat("\n Potential outliers(Local) = ", Loutlier.potential)
  
  if( is.na(Loutlier.potential[1]) == FALSE )  
  {
    originGroup[Loutlier.potential] = 0
    origin.cluster = originGroup[originGroup > 0]
    cluster.id = clusterGroup[originGroup > 0]
    km.data = kdata[-Loutlier.potential,c(selected.vars)]
    adjRand.out = adjwithout.outlier(km.data, cls.num, origin.cluster,cluster.id)
    # Show the Ratio of betweenSS/totSS without outlier
    ratio.val = ratio.calculate(km.data, cluster.id)
  }
  else
  {
    adjRand.out = adjwithout.outlier(kdata[,c(selected.vars)], cls.num, originGroup, clusterGroup)
    # Show the Ratio of betweenSS/totSS without outlier
    ratio.val = ratio.calculate(kdata[,c(selected.vars)], clusterGroup)
  }
  
  cat("\n Original SSB/SST = ", round(ratio.val,4), "\n")
  cat("\n Adjusted Rand Index = ", adjRand.out$adjRand,"\n")
  cat("\n Confuion Matrix \n") 
  print(adjRand.out$clsTable)
  
}


adjwithout.outlier <- function(km.data, clsGroup, origin.cluster, cluster.id)
{
  #----------------------------------------------------------#
  # K-means with selected variables and without outliers     #
  #----------------------------------------------------------#
  
  Gmean <- matrix(0, clsGroup, ncol(km.data))
  
  clsInfo = cluster.id
  clsdata <- cbind(clsInfo, km.data)
  
  
  Glength <- array(0, clsGroup)
  
  for(i in 1:clsGroup) { Glength[i] <- length(clsInfo[clsInfo == i]) }
  
  for(i in 1:clsGroup)
  {    gdata <- clsdata[clsdata[,1] == i,]
  if( Glength[i] == 1 )
  {  Gmean[i,] <- gdata[-1] }
  else
  { 
    Gm <- apply(gdata, 2, 'mean')
    Gmean[i,] <- Gm[-1] 
  }
  }
  
  km.out =  kmeans(km.data, Gmean)
  
  cls.id = km.out$cluster
  cls.size = km.out$size
  origin.length = length(origin.cluster)
  cluster.length = length(cls.id)
  
  adjRand2 = adjustedRandIndex(cls.id, origin.cluster) 
  adjRand2 = round(adjRand2,4)
  
  cls.table = table(origin.cluster, cls.id)
  adjRand.out = list(adjRand = adjRand2, clsTable = cls.table)
  return(adjRand.out)
}

pkgInstall.check <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

SVOKmeans <- function(data=datafile)
{
  pkgInstall.check("mclust")
  pkgInstall.check("sampling")
  pkgInstall.check("chemometrics")
  
  library(mclust)
  library(sampling)  
  library(chemometrics)
  
  no.col <- ncol(data)
  no.row <- nrow(data)
  
  
  
  # ------------------------------------------------------------#
  # Step 0-1 : Data Transformation                              #
  # ------------------------------------------------------------#
  cat("\n -----------------------------------------")
  cat("\n Step 0-1: Standardize Variables ?       ")
  cat("\n   1. 0-1 Transform  2. Z-Score  3. Raw Data")
  cat("\n -----------------------------------------")
  cat("\n   Select(default:1) : ")
  s.type <- readline()
  if(s.type == "" ) s.type <- 1
  s.type <- as.integer(s.type)
  
  if( s.type == 1)
  {   maxX = apply(data,2, max)
  minX = apply(data,2, min) 
  kdata <- scale(data, center=minX, scale=maxX-minX)
  }
  else if(s.type == 2) 
    kdata <- scale(data, center=TRUE, scale=TRUE)
  else
    kdata <- data 
  
  kdata = data.frame(kdata)
  
  # ------------------------------------------------------------#
  # Step 0-2 : Sampling                                         #
  # ------------------------------------------------------------#
  cat("\n -----------------------------------------------------")
  cat("\n Step 0-2: 1. Sampling 2. Full Data : # of data=", no.row," ") 
  cat("\n -----------------------------------------------------")
  cat("\n   Select(default=1): ")
  data.select <- readline()
  if(data.select == "") data.select <- 1
  data.select <- as.integer(data.select)
  if(data.select != 2) data.select <- 1
  
  if(data.select == 1 ) 
    sam.data = sample.data(kdata)
  else 
    sam.data = kdata
  
  # ------------------------------------------------------------#
  # Step 0-3 : Parameter for Deciding the number of clusters    #
  # ------------------------------------------------------------#
  cat("\n------------------------------------------------------------------------")
  cat("\n Step 0-3: Mojena's k for deciding the number of clusters(def=1.25): ")
  wk = readline()
  if(wk == "") wk=1.25 
  else 
    wk = as.numeric(wk)
  
  cat("\n------------------------------------------------------------------------")
  cat("\n Step 0-4: Parameter for Outlier                   ") 
  cat("\n       - Mahal. distance(1=default) Robust Mahal. distance(2) : ")
  mahal.select = readline()
  if(mahal.select == "") mahal.select = 1  
  else 
    mahal.select = as.numeric(mahal.select)
  
  if(mahal.select == 2) 
  {  cat("         Robust Function : cov.mve(1=def) cov.rob(2)  cov.mcd(3) covMcd(4) : ")
    rtype.sel = readline()
    if(rtype.sel == "") rtype.sel=1 
    else
      rtype.sel = as.numeric(rtype.sel)
  }
  else
    rtype.sel = 1
  
  cat("       - Chisq-quantile for outlier(def.=0.975): ")
  q.value = readline()
  if(q.value == "") q.value = 0.975  
  else 
    q.value = as.numeric(q.value)
  
  cat("       - Ratio of cluster small size for outlier(def=0.03): ")
  cls.ratio = readline()
  if(cls.ratio == "") cls.ratio = 0.03  
  else 
    cls.ratio = as.numeric(cls.ratio)
  
  cat("------------------------------------------------------------------------\n")
  
  # -----------------------------------------------------------------
  # Step 1 : Automated K-means for each variables                  
  # -----------------------------------------------------------------
  # Initialize 
  no.col <- ncol(kdata)
  id.label = rownames(kdata)  # identification number of observations
  kno.row <- nrow(kdata)
  ino.row <- nrow(sam.data)
  Pj = matrix(0, ncol=no.col, nrow=kno.row)
  y  = array(0, kno.row)
  PjGroup = array(0, no.col)
  
  for(j in 1:no.col)
  { 
    e.data = sam.data[,j]
    InitRes <- WardClustering(e.data,wk)
    clsGroup <- InitRes$clsGroup
    PjGroup[j] <- clsGroup
    initCenter <- InitRes$clsCenter
    
    kmeans.v <- kmeans(kdata[,j], initCenter)
    Pj[,j] <- kmeans.v$cluster           
  }
  
  cat("\n < Variable Selection and Outlier Detection for K-means > \n")
  
  cat("\n Step 1 : Size of Cluster for each variable \n")
  print(PjGroup)
  
  #-------------------CHECK ----------------------------------------- 
  wm.a = which.max(table(PjGroup)) 
  wm.a = names(wm.a)
  wm.a = as.integer(wm.a)
  initKmeansGroup = wm.a 
  
  # cat("        : Initial K-means Group = ", initKmeansGroup )
  #------------------------------------------------------------------
  
  # ------------------------------------------------------------------------------
  # Step 2 : Adjusted Rand Index for each pair of variables and find highest top 5                 
  # ------------------------------------------------------------------------------
  
  adj.mat = matrix(0, ncol=no.col, nrow=no.col)
  
  no.d = no.col - 1
  for( j in 1:no.d )
    for(k in (j+1):no.col)
      adj.mat[j,k] = adj.mat[k,j] = adjustedRandIndex(Pj[,j], Pj[,k]) 
  
  # cat("\n ===  Adjusted Rand Index ===\n")
  # print(adj.mat)
  
  if( no.col == 2 ) adj.row = 1
  else if (no.col == 3) adj.row = 3
  else adj.row=5
  
  top5.adj = matrix(0, ncol=2, nrow=adj.row)
  top5.adj.index = array(0,adj.row)
  
  for(i in 1:adj.row)
  { 
    max.number = which.max(adj.mat)
    max.v1 = max.number %% no.col
    max.v2 = max.number %/% no.col  
    
    if(max.v1 == 0 ) 
      max.v1 = no.col
    else
      max.v2 = max.v2 +1
    
    # cat("\n num max.v1 and v2 = (", max.number, ",", max.v1, ",", max.v2, ")")
    top5.adj[i,1] = max.v2
    top5.adj[i,2] = max.v1
    
    top5.adj.index[i] = adj.mat[max.v1, max.v2]
    adj.mat[max.v1, max.v2] = adj.mat[max.v2, max.v1] = -30
    
  }
  
  top5.adj.index = top5.adj.index[top5.adj.index > 0]
  adj.row = length(top5.adj.index)
  top5.adj = top5.adj[c(1:adj.row),]  
  top5.ssb.sst = array(0, adj.row)
  
  cat("\n\n Step 2 : Highest adjusted Rand Index(Up to Top 5) \n")
  print(top5.adj)
  
  #---------------------------------------------------------------------------------
  # Step 3 : Calculate ratio of SSB/SST
  #---------------------------------------------------------------------------------
  
  for(i in 1:adj.row)
  {
    vars = top5.adj[i,]
    # cat("\n vars = ", vars)
    e.data = sam.data[,c(vars)]
    k.data = kdata[,c(vars)]
    kmeans.out = TwoStepKmeans(e.data, k.data,wk)
    ss.ratio = kmeans.out$betweenss / kmeans.out$totss
    top5.ssb.sst[i] = ss.ratio      
  }
  
  cat("\n Step3 : Ratio of SSB/SST (Up to Top 5) \n")
  print(round(top5.ssb.sst,4))
  
  
  #---------------------------------------------------------------------------------
  # Step 4 : Select two variables which have the highest ratio of SSB/SST
  #---------------------------------------------------------------------------------
  
  top.index = which.max(top5.ssb.sst)
  top.vars = top5.adj[top.index,]
  all.vars = c(1:no.col)
  selected.vars = top.vars
  unselected.vars = all.vars[-selected.vars]
  
  cat("\n Step4 : First Selected Var's = (", selected.vars,")" )
  
  
  #---------------------------------------------------------------------------------
  # Step 5 : Run K-means using the selected variables in Step 4
  #---------------------------------------------------------------------------------
  
  pselected.vars = sort(selected.vars)
  e.data = sam.data[,c(pselected.vars)]
  k.data = kdata[,c(pselected.vars)]
  kmeans.out = TwoStepKmeans(e.data, k.data,wk)
  
  cat("\n\n Step 5 : Results of first selected var's ")
  
  # Outlier Detection Procedure with selected two variables
  cls.id = kmeans.out$cluster
  cls.size = kmeans.out$size
  cls.num = length(kmeans.out$size)
  cls.center = kmeans.out$centers
  
  cat("\n        Selected Var's = (", selected.vars,")" )
  unselected.num = length(unselected.vars)
  if(unselected.num > 0 ) cat("\n        UnSelected Vars = (", unselected.vars, ")" )
  else cat("\n        UnSelected Vars = ( - )" )
  cat("\n        Number of Cluster = ", cls.num)
  cat("\n        Cluster Sizes = ", cls.size)
  
  outvar2 = Outlier.detection.fn(k.data, cls.center, cls.id, cls.num, id.label,q.value, cls.ratio, mahal.select, rtype.sel) 
  
  Loutlier = outvar2$Loutlier
  Ldistance = outvar2$Ldistance
  Ldtype = outvar2$Ldtype
  Goutlier = outvar2$Goutlier
  Gdistance = outvar2$Gdistance
  Gdtype = outvar2$Gdtype
  ss.ratio =  kmeans.out$betweenss / kmeans.out$totss
  cutoff = outvar2$cutoff
  
  PrintOutlier.ftn(Loutlier, Ldistance, Ldtype,Goutlier, Gdistance, Gdtype, cutoff, ss.ratio, mahal.select) 
  
  YPj = kmeans.out$cluster
  
  #---------------------------------------------------------------------------------
  # Step 6 :  kmeans for unselected variables and  Compute Adjusted Rand Index 
  #           between Selected Vars and Unselected Vars 
  #---------------------------------------------------------------------------------
  
  y.u.adj = array(0, no.col)
  
  for(j in c(unselected.vars) )
  {
    num.unique = length(unique(kdata[,j]))
    if( num.unique < cls.num ) Ugroup = kdata[,j]
    else
    { kmeans.v <- kmeans(kdata[,j], cls.num)
    Ugroup <- kmeans.v$cluster           
    } 
    y.u.adj[j] = adjustedRandIndex(YPj, Ugroup)  
  }
  
  cat("\n\n Step 6 : Adjusted Rand Index between Y and Unselected var's \n")
  print(round(y.u.adj,4))
  
  u.adj.max = max(y.u.adj)
  u.adj.which = which.max(y.u.adj)
  
  
  #---------------------------------------------------------------------------------
  # Step 7 :  If max rand index > critical value, add selected variable, or Stop 
  #---------------------------------------------------------------------------------
  
  cat("\n Step 7 : Repeating procedure for adding var's ")
  critical.value = 0.05 
  eta = u.adj.max
  Gfac = 0.5 
  u.leng = length(unselected.vars) 
  
  
 tryCatch({ while(u.adj.max >= critical.value | u.adj.max >= eta*Gfac && u.leng > 0 ) 
  { 
    selected.vars = all.vars[c(selected.vars, u.adj.which)] 
    unselected.vars = all.vars[-selected.vars] 
    eta = u.adj.max
    u.leng = length(unselected.vars) 
    kresult = RunStep5.Step6.rtn(sam.data, kdata, selected.vars, unselected.vars,wk, id.label, q.value, cls.ratio, mahal.select,rtype.sel)
    u.adj.max = kresult$adjmax
    u.adj.which = kresult$adjwhich
    cls.id = kresult$cls.id
    cls.size = kresult$cls.size
    Loutlier.potential = kresult$Loutlier
    Loutlier.distance = kresult$Ldistance
    Ldtype = kresult$Ldtype
    Goutlier.potential = kresult$Goutlier
    Goutlier.distance = kresult$Gdistance
    Gdtype = kresult$Gdtype
    cutoff = kresult$cutoff
    ss.ratio = kresult$ratio
    # cat("\n adj.max = ", round(u.adj.max,4), "which= ", u.adj.which)
  }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  #---------------------------------------------------------------------------------
  # Step 8 :  Last Kmeans Results using selected variables  
  #---------------------------------------------------------------------------------
  
  cat("\n\n Step 8 : Last K-Means Results Using Selected Variables ")
  cat("\n        Selected Vars = (", selected.vars,")" )
  
  unselected.num = length(unselected.vars)
  if(unselected.num > 0 ) cat("\n        UnSelected Vars = (", unselected.vars, ")" )
  else cat("\n        UnSelected Vars = ( - )" )
  
  cat("\n        Number of Cluster = ", length(cls.size) )
  cat("\n        Cluster Size = ", cls.size)  
  
  PrintOutlier.ftn(Loutlier.potential, Loutlier.distance, Ldtype, Goutlier.potential, Goutlier.distance, Gdtype, cutoff, ss.ratio, mahal.select) 
  
  uniq.CL = unique(Ldtype)
  uniq.GL = unique(Gdtype)
  
  if(any(uniq.CL>2) || any(uniq.GL>2) )
  {
    cat("\n ----------------------------------------------------------------------------------")
    cat("\n (*)Comment in 'Outlier Type'  ")
    cat("\n    (1) sn :  small size clusters(ratio=",cls.ratio,")")
    cat("\n    (2) sp :  number of cases in cluster is too small compared to variables  ") 
    cat("\n    (3) md :  Mahalanobis distance is provided in cases Robust Mahal. distcance ") 
    cat("\n              can not be calculated                                             ")
    cat("\n    (4) -  :  (Robust) Mahalanobis distance rquired                             ")
    cat("\n ----------------------------------------------------------------------------------")
  }
  
  cat(" -----------------------------------------------------------------------------------")
  cat("\n   Adjusted Rand Index ? - Simulated Data(1), Real Data(2), None(3=default) : ")
  AduCal = readline()
  if(AduCal == "") AduCal=3
  AduCal = as.integer(AduCal)
  cat(" -----------------------------------------------------------------------------------")
  
  if (AduCal == 1) AdjustedSim.Calculate(kdata, id.label, cls.size, cls.id, all.vars, selected.vars,unselected.vars, Loutlier.potential, Goutlier.potential ) 
  else if (AduCal == 2) AdjustedRdata.Calculate(kdata, cls.size, cls.id, selected.vars, Loutlier.potential, Goutlier.potential) 
  
  cat("\n ===  End of Processing  ===\n")
  kmeans.out[['selected.vars']] <- selected.vars
  kmeans.out[['names']] <- names(kdata)
  return(kmeans.out)
}
