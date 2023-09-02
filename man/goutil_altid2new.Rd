\name{goutil$altid2new} 
\alias{goutil_altid2new} 
\alias{goutil$altid2new} 
\title{Convert old alternative GO ids to their new counterpart} 
\description{ 
This function returns the names for a vector of given GO ids.  
Before the function can work first a GO obo file must be read in using the function `goutil$new`. 
} 
\usage{ goutil_altid2new(goid) } 
 
\arguments{ 
  \item{goid}{vector of Gene Ontology ids} 
} 
\value{vector of names belonging to the given GO id(s)} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
goutil$altid2new(c("GO:0000001","GO:0061637")) # first should stay unchanged 
} 
 
