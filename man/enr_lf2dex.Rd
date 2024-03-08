\name{enr$lf2dex} 
\alias{enr_lf2dex} 
\alias{enr$lf2dex} 
\title{Create T/F table for a certain log-fold difference} 
\description{ 
  Create a TRUE/FALSE data frame to indicate if the given two log-fold expressions 
  against a control for two different conditions are differentially expressed 
} 
\usage{ enr_lf2dex(x, threshold.control=1, threshold.groups=0.5) } 
 
\arguments{ 
 \item{x}{a two column data frame with log-fold change values for two groups against a control} 
 \item{threshold.control}{the numerical threshold for the comparison against the control},  
 \item{threshold.groups}{the numerical threshold for the comparison between the groups}against the control,  
} 
\details{ 
  Usually such an analysis can be down using a simple check against a given 
  log-fold change against a control for instance, if log2 values are compared 
  data>1, however for the difference set between the two groups, it can happen 
  that one group has a log-fold change against the control, let's say of 1.3 
  and the the other group a value of 0.9 against the control. So between each 
  other they are not differentially expressed. The TRUE/FALSE value in this case 
  for a  simple check against > 1 would be TRUE (G1) against the control and  
  FALSE (G2) against the control. 
} 
\value{data frame with two columns with tTRUE or FALSE values for both comparisons} 
\examples{ 
sdata=data.frame( 
G1=c(1.2,1.2,1.6,1.2,0.8,-1.6,-1.2,-1.2,-0.1,-0.3,-0.8,0.4,-3.6,-2.6), 
G2=c(0.9,0.1,1.6,1.8,2.5,-1.4,-0.1,-0.3,-1.2,-1.2,-1.7,1.7,-1.4,-1.4)) 
rownames(sdata)=LETTERS[1:14] 
sdata 
## up-regulated 
sdata > 1 
enr$lf2dex(sdata,threshold.control =  1) 
## down-regulated 
sdata < -1 
enr$lf2dex(sdata,threshold.control = -1) 
dex=enr$lf2dex(sdata,threshold.control =  1) 
table(dex[,1],dex[,2]) 
enr$lf2dex(sdata,threshold.control =  1, threshold.groups=0.3) 
} 
 
