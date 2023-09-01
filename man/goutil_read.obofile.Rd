\name{goutil$read.obofile} 
\alias{goutil_read.obofile} 
\alias{goutil$read.obofile} 
\title{ Read the given Gene Ontology Obo file } 
\description{ 
  This function does the parsing of the obofile, should be used mainly indirectly via `goutil$new`. 
} 
\usage{ goutil_read.obofile(obofile) } 
 
\arguments{ 
  \item{obofile}{ A GO obo file } 
} 
\details{ 
  This function read the actual obofile. 
} 
\value{list with the following components: 
\itemize{ 
 \item names - with columns id, nsp, and name 
 \item slims - with columns id, slim 
 \item obsoletes - with ids which are obsolete 
 \item tree - with columns child, parent, relation 
} 
} 
\examples{ 
res=goutil$read.obofile("2023-01-01-go.obo") 
sort(names(res)) 
} 
 
