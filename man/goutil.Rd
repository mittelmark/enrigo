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
\item \code{\link[enrigo:goutil_new]{goutil$new}}{initialize the required GO data using the given obo file} 
\item \code{\link[enrigo:goutil_altid2new]{goutil$altid2new}}{convert old alternative GO ids to their new counterpart} 
\item \code{\link[enrigo:goutil_getChildren]{goutil$getChildren}}{get the child nodes of a given GO id} 
\item \code{\link[enrigo:goutil_getEntry]{goutil$getEntry}}{get the GO entry in standard text for a given GO id} 
\item \code{\link[enrigo:goutil_getName]{goutil$getName}}{get the name of (a) given GO ids} 
\item \code{\link[enrigo:goutil_getNamespace]{goutil$getNamespace}}{get the namespace of (a) given GO ids} 
\item \code{\link[enrigo:goutil_getOboVersion]{goutil$getOboVersion}}{get the actual version date for the current obo file} 
\item \code{\link[enrigo:goutil_getParent]{goutil$getParent}}{get the parent node(s) for a given GO id} 
\item \code{\link[enrigo:goutil_getSlims]{goutil$getSlims}}{get GO-slims, or the GO ids for a given slim} 
\item \code{\link[enrigo:goutil_getStats]{goutil$getStats}}{statistics for the three main namespaces and the number of obsolete terms} 
\item \code{\link[enrigo:goutil_getTree]{goutil$getTree}}{reads recursively all parents nodes for a given GO id} 
\item \code{\link[enrigo:goutil_getTreeMatrix]{goutil$getTreeMatrix}}{adjacency matrix for the given ids in tree} 
\item \code{\link[enrigo:goutil_isChild]{goutil$isChild}}{check if a given GO id is a direct child of a given parent id} 
\item \code{\link[enrigo:goutil_isParent]{goutil$isParent}}{check if a given GO id is a direct parent of a given child id} 
\item \code{\link[enrigo:goutil_kroki]{goutil$kroki}}{create a GO tree graph using the kroki webservie} 
\item \code{\link[enrigo:goutil_quickGO]{goutil$quickGO}}{create a GO tree graph using the EBI QuickGO webservie} 
\item \code{\link[enrigo:goutil_read.obofile]{goutil$read.obofile}}{ reads the given obofile and returns the result} 
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
