\name{goutil$getParent} 
\alias{goutil_getParent} 
\alias{goutil$getParent} 
\title{Return the parent GO id(s) for given GO id(s)} 
\description{ 
This function returns a vector of parent GO id(s) for all given GO id(s). 
Before the function can work first a GO obo file must be read in using the function `goutil$new`. 
} 
\usage{ goutil_getParent(goid,type='all') } 
 
\arguments{ 
  \item{goid}{single Gene Ontology ids} 
  \item{type}{which type of parents should be returned, either 'all', 'is_a' or 'part_of', default: 'all'} 
} 
\details{ 
 This function returns the parent GO ids for a vector of given GO ids. This function can be used 
 for explore the tree of all parent ids for a given GO id. 
 Before the function is called you should initialize the goutil environment with the  
 function `goutil$new(obofile)` 
} 
\value{vector of parent ids belonging to the given GO id} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
goutil$getParent("GO:0000001") 
} 
 
