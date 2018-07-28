.findPos=function(target,sc)
{
  outcome=NULL;
  if(length(target)>0)
  {
    outcome=integer(length(target));
    for(i in 1:length(outcome))
      outcome[i]=which(sc==target[i]);
  }
  return(outcome);
}

.fbplotDetect<-function(data,depth,shapeOutlier)
{
  nCurve=dim(data)[2];
  nCentral=ceiling(nCurve/2);

  remain=1:nCurve;
  if(!is.null(shapeOutlier))
  {
    data=data[,-shapeOutlier]
    depth=depth[-shapeOutlier]
    remain=remain[-shapeOutlier];
  }
  index=order(depth,decreasing=T);

  center=data[,index[1:nCentral]];
  inf=apply(center,1,min);
  sup=apply(center,1,max)
  dist=1.5*(sup-inf);
  upper=sup+dist;
  lower=inf-dist;

  outly=(data<=lower)+(data>=upper);
  outpoint=which(colSums(outly)>0);
  return(remain[outpoint]);
}

.boxDetect<-function(depth,boxFactor)
{
  result=boxplot(depth,range=boxFactor,plot=F);
  out=result$out[result$out<mean(depth)];
  return(.findPos(out,depth));
}

#' Return TVD and MSS for functional data
#'
#' @param nCurve A scalar, the number of curves
#' @param nPoint A scalar, the number of discrete time points
#' @param data A matrix with dimension nCurve by nPoint, the functional data
#' @return TVD and MSS
TVDMSS=function(data,nCurve,nPoint)
{
  if(nCurve<3) stop("the number of curves must be greater than 2")
  pointwiseRank=t(apply(data,1,rank));

  shapeVar=matrix(0,nPoint,nCurve);
  totalVar=pointwiseRank*(nCurve-pointwiseRank)/nCurve/nCurve;

  pointwiseMedian=apply(data,1,median);
  for(i in 2:nPoint)
  {
    medianAtI=pointwiseMedian[i];
    for(j in 1:nCurve)
    {
      dataAtI=data[i,];shift=dataAtI[j]-medianAtI;
      dataAtI_=data[i-1,];
      dataAtI[j]=medianAtI;dataAtI_[j]=dataAtI_[j]-shift;

      belowAtI_=dataAtI_<=dataAtI_[j];
      belowAtI=dataAtI<=medianAtI;

      below_=sum(belowAtI_);
      below_below=sum(belowAtI_ & belowAtI);
      above_below=sum((!belowAtI_)&belowAtI);

      part1=below_below^2/below_;
      if(below_==nCurve) part2=0 else part2=above_below^2/(nCurve-below_);

      shapeVar[i,j]=(part1+part2)/nCurve-(nCurve/2)^2/nCurve/nCurve;
    }
  }

  v=abs(apply(data,2,diff))
  v_=colSums(v);
  v=t(t(v)/v_);

  TVD=colMeans(totalVar);

  MSS=shapeVar[2:nPoint,]/0.25;
  MSS=colSums(MSS*v);

  return(list(TVD=TVD,MSS=MSS));
}

detectOutlier<-function(data,nCurve,nPoint,empFactor)
{
  if(nCurve<3) stop("the number of curves must be greater than 2")
  # compute TVD and MSS
  myDepth=TVDMSS(data,nCurve,nPoint);
  TVD=myDepth$TVD;
  MSS=myDepth$MSS;

  # find outliers
  magOutlier=NULL;shapeOutlier=NULL;
  shapeOutlier=.boxDetect(MSS,empFactor)
  magOutlier=.fbplotDetect(data,TVD,shapeOutlier);

  return(list(outlier=sort(c(shapeOutlier,magOutlier)),sOut=shapeOutlier,mOut=magOutlier,TVD=TVD,MSS=MSS));
}
