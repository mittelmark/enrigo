\name{uniprot$ensg2id} 
\alias{uniprot_ensg2id} 
\alias{uniprot$ensg2id} 
\title{Extract mappings between uniprot identifier and Ensembl gene indentifier } 
\description{ 
  This function creates a mapping file between Uniprot identifiers and Ensembl gene identifiers. 
} 
\usage{ uniprot_ensg2id(infile,outfile) } 
 
\arguments{ 
  \item{infile}{UniProt file downloaded with `uniprot$download`} 
  \item{outfile}{file where to store the data} 
} 
\details{ 
  This function needs up to several minutes in dependence from the file size. 
} 
\value{NULL} 
\examples{ 
# get all species 
file=uniprot$download("mammals", 
 folder=file.path(path.expand("~"),"data")) 
id2ensgfile=file.path(path.expand("~"),"data","id2ensg.tab") 
uniprot$id2ensg(file,id2ensgfile) 
head(read.table(id2ensgfile,sep="\t",quote="")) 
} 
 
