\name{enr$new} 
\alias{enr_new} 
\alias{enr$new} 
\title{Initialize the annotation data } 
\description{ 
  This function reads in the Gene Ontology annotation file.  
   
} 
\usage{ enr_new(gaffile=NULL) } 
 
\arguments{ 
  \item{gaffile}{GO annotation filename which should be used to perform the  
   enrichment analysis, if not given the last downloaded file, using the function `enr$gaf` is used, default: NULL } 
} 
\details{ 
  This function allows you to initialize the annotation data which contain the mapping between 
  the GO ids and the gene ids for a given species. 
} 
\value{NULL} 
\examples{ 
# get all species 
enr$gaf("Apis mellifera") 
enr$new() 
} 
 
