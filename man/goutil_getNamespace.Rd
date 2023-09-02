\name{goutil$getNamespace} 
\alias{goutil_getNamespace} 
\alias{goutil$getNamespace} 
\title{Return the namespaces(s) for given GO id(s)} 
\description{ 
This function returns a vector of namespaces for all given GO id(s). 
Before the function can work first a GO obo file must be read in using the function `goutil$new`. 
} 
\usage{ goutil_getNamespace(goid) } 
 
\arguments{ 
  \item{goid}{vector of Gene Ontology id(s)} 
} 
\details{ 
 This function returns the namespaces with single letter codes for a vector of given GO ids. 
 Before the function is called you should initialize the goutil environment with the  
 function `goutil$new(obofile)` 
} 
\value{vector of namespace abbreviations belonging to the given GO id(s)} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
goutil$getNamespace(c("GO:0003675","GO:0008150")) 
} 
 
