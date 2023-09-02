\name{goutil$kroki} 
\alias{goutil_kroki} 
\alias{goutil$kroki} 
\title{Return an URL for a PlantUML graph done using the kroki webservice} 
\description{ 
This function can be used to embed images for Gene Ontology trees into an HTML page. 
} 
\usage{ goutil_kroki(goid,plot=FALSE) } 
 
\arguments{ 
  \item{goid}{single Gene Ontology id} 
  \item{plot}{should the image plotted into a standard R plot, requires the png package, default: FALSE} 
} 
\value{returns a URL which can be embedded into a Markdown URL for instance} 
\examples{ 
goutil$kroki("GO:0003676") 
} 
 
