\name{goutil}  
\alias{goutil} 
\alias{goutil-class} 
\title{ Environment object with functions to work with Gene Ontology terms } 
\description{ 
This is a collection of useful functions to work with Gene Ontology terms and 
annotations. Enrichment analysis for gene lists and graphical visualizations are as well provided.#' Below you find the documentation for the following functions  
} 
\section{Methods}{ 
\itemize{ 
\item \code{\link[enrigo:goutil_new]{goutil$new(obofile)}}{initialize the required GO data using the given obo file} 
\item \code{\link[enrigo:goutil_read.obofile]{goutil$read.obofile(obofile)}}{ reads the given obofile and returns the result} 
} 
} 
\examples{ 
ls(goutil) 
# get all species 
fname=goutil$download(2023) 
print(fname) 
res=goutil$new(fname) 
names(res) 
}  
