\name{goutil$getEntry} 
\alias{goutil_getEntry} 
\alias{goutil$getEntry} 
\title{Return the the GO entry as standard text for a given GO id} 
\description{ 
This function returns the GO entry as a character vector with one line per element. 
Before the function can work first a GO obo file must be read in using the function `goutil$new`. 
} 
\usage{ goutil_getEntry(goid) } 
 
\arguments{ 
  \item{goid}{single Gene Ontology ids} 
} 
\value{text entry belonging to the given GO id} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
writeLines(unlist(lapply(goutil$getEntry("GO:0000001"),strwrap, 80, exdent=6))) 
} 
 
