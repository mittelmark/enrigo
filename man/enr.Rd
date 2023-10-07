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
\item \code{\link[enrigo:enr_annotation]{enr$annotation(genes)}}{return GO annotations for the given gene identifiers} 
\item \code{\link[enrigo:enr_enrichment]{enr$enrichment(fullset,subset)}}{perform an enrichment analysis} 
\item \code{\link[enrigo:enr_gaf]{enr$gaf()}}{query, download or initialize a Gene Ontology association file} 
\item \code{\link[enrigo:enr_download]{enr$download(species)}}{dowmload an annotation file for the given species} 
\item \code{\link[enrigo:enr_new]{enr$new(filename)}}{initialize the GO annotation} 
\item \code{\link[enrigo:enr_name2id]{enr$name2id(names)}}{convert free text terms into identifiers} 
\item \code{\link[enrigo:enr_symbol2loc]{enr$symbol2loc(symbol)}}{convert gene esymbols into LOC identifiers} 
} 
} 
\examples{ 
# get all species 
res=enr$gaf() 
sort(names(res)) 
}  
