\name{goutil$getSlims} 
\alias{goutil_getSlims} 
\alias{goutil$getSlims} 
\title{Return all GO-slims, or the GO ids for a given slim and possibly within the given namespace} 
\description{ 
This function returns either all existing GO slim names or if 
if a slim name is given all GO ids, belonging to this slim, if a namespace is given as well the GO ids are checked for the 
given namespace as well. 
} 
\usage{ goutil_getSlims(slim,nsp) } 
 
\arguments{ 
  \item{slim}{optional slim name, default: NULL} 
  \item{nsp}{optional namespace like f (molecular function), c (cellular component), p (biological process), if slim is given and the nsp only GO ids from this namespace will be returned, default: NULL} 
} 
\value{vector of slim names or if a slim name is given all GO ids belonging to this slim} 
\examples{ 
goutil$getSlims() 
goutil$getSlims(slim="drosophila",nsp="c") 
for (n in c('c','f','p')) { 
  print(paste(n,"=",length( goutil$getSlims("drosophila",nsp=n)))) 
} 
} 
 
