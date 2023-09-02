\name{goutil$getChildren} 
\alias{goutil_getChildren} 
\alias{goutil$getChildren} 
\title{Return the child nodes of the given GO id} 
\description{ 
 Return the GO id(s) which have the given go id as a 'is_a' or 'part_of' parent. 
 Before the function can work first a GO obo file must be read in using the function `goutil$new`. 
} 
\usage{ goutil_getChildren(goid) } 
 
\arguments{ 
  \item{goid}{single Gene Ontology ids} 
} 
\value{vector of child ids having the given GO id as parent} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
goutil$getChildren("GO:0003674") 
} 
 
