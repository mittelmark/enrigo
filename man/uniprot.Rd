\name{uniprot}  
\alias{uniprot} 
\alias{uniprot-class} 
\title{ Environment object with functions to prepare UniProt data for enrichment analysis } 
\description{ 
The functions within the uniprot environment allow to extract data from UniProt files to do Gene Ontology term enrichment analysis for a given list 
of genes against a full list of genes from the same organism. 
} 
\section{Methods}{ 
\itemize{ 
\item \code{\link[enrigo:uniprot_download]{uniprot$download("taxdivision")}}{download a SwissProt file from the UniProt database from the taxonomic division section} 
\item \code{\link[enrigo:uniprot_id2go]{uniprot$id2go(infile,outfile)}}{extract mappings between UniProt identifier and GO identifiers} 
\item \code{\link[enrigo:uniprot_ensg2id]{uniprot$ensg2id2(infile,outfile)}}{extract mappings between UniProt identifiers and Ensembl identifiers} 
\item \code{\link[enrigo:uniprot_mapid2go]{uniprot$mapid2go(upid2gofile,idfile,outfile,evidence="ALL")}}{create mappings between identifiers and GO identifiers using UniProt identifiers} 
} 
} 
\examples{ 
# download mammmals SwisProt database file 
file=uniprot$download("human", 
  folder=file.path(path.expand("~"),"data")) 
file 
}  
