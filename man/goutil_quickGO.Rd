\name{goutil$quickGO} 
\alias{goutil_quickGO} 
\alias{goutil$quickGO} 
\title{Return an URL for a QuickGO graph using the EBI webservice} 
\description{ 
This function can be used to embed images for Gene Ontology trees into a Markdown or HTML document. 
} 
\usage{ goutil_quickGO(goid,key=FALSE,children=FALSE,plot=FALSE) } 
 
\arguments{ 
  \item{goid}{single Gene Ontology id} 
  \item{key}{show the legend, default: FALSE} 
  \item{children}{show the child terms of a GO id} 
  \item{plot}{should the image plotted into a standard R plot, requires the png package, default: FALSE} 
} 
\value{returns a URL which can be embedded into a Markdown URL for instance} 
\examples{ 
goutil$quickGO("GO:0000001") 
goutil$quickGO("GO:0003674",children=TRUE) 
} 
 
