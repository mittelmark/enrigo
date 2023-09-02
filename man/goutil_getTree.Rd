\name{goutil$getTree} 
\alias{goutil_getTree} 
\alias{goutil$getTree} 
\title{Return all parent GO id(s) recursively for given GO id(s)} 
\description{ 
This function returns a vector of all GO ids of parents and their parents 
  an so on, so recursively. 
} 
\usage{ goutil_getTree(goid) } 
 
\arguments{ 
  \item{goid}{single Gene Ontology ids} 
} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
  for (id in goutil$getTree("GO:0003676")) {  
     cat(id,"\t",goutil$getName(id),"\n")  
  } 
} 
 
