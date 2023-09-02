\name{goutil$getName} 
\alias{goutil_getName} 
\alias{goutil$getName} 
\title{Return the name(s) for given GO id(s)} 
\description{ 
This function returns a vector of names for all given GO id(s). 
Before the function can work first a GO obo file must be read in using the function `goutil$new`. 
} 
\usage{ goutil_getName(goid) } 
 
\arguments{ 
  \item{goid}{vector of Gene Ontology ids} 
} 
\details{ 
  This function returns the names for a vector of given GO ids. Before the function 
  is called you should initialize the goutil environment with the function `goutil$new(obofile)` 
} 
\value{vector of names belonging to the given GO id(s)} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
goutil$getName(c("GO:0000001","GO:0000002")) 
goutil$getName("GO:0003675") 
} 
 
