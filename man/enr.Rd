\name{enr}  
\alias{enr} 
\alias{enr-class} 
\title{ Environment object with functions do do enrichment analysis for GO terms } 
\description{ 
The functions within the enr environment allow to do Gene Ontology term enrichment analysis for a given list 
of genes against a full list of genes from the same organism. 
} 
\section{Methods}{ 
\itemize{ 
\item \code{\link[enrigo:enr_gaf]{enr$gaf(filename)}}{query, download or initialize a Gene Ontology association file} 
} 
} 
\examples{ 
# get all species 
res=enr$gaf() 
sort(names(res)) 
}  
