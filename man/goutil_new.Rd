\name{goutil$new} 
\alias{goutil_new} 
\alias{goutil$new} 
\title{ Initializes the goutil object with a new obofile } 
\description{ 
  This function should be used first to initialize the internal data tables and information. 
} 
\usage{ goutil_new(obofile,cache=TRUE) } 
 
\arguments{ 
  \item{obofile}{GO obo filename } 
  \item{cache}{should the results of the parins being cached, default: TRUE } 
} 
\details{ 
  This function allows you to intialize your annotation file which contains the mapping between 
  the GO ids and the gene ids for a given species. For getting a list of available species and for the  
  download a internet connection is required. 
} 
\value{list with the folloowing components: 
\itemize{ 
 \item names - with columns id, nsp, and name 
 \item slims - with columns id, slim 
 \item obsoletes - with ids which are obsolete 
 \item tree - with columns child, parent, relation 
} 
} 
\examples{ 
res=goutil$new("2023-01-01-go.obo") 
sort(names(res)) 
} 
 
