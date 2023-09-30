\name{uniprot$id2go} 
\alias{uniprot_id2go} 
\alias{uniprot$id2go} 
\title{ Extract mappings between uniprot identifier and GO identifiers } 
\description{ 
  This function creates a mapping file between Uniprot identifiers and GO terms. 
} 
\usage{ uniprot_id2go(infile,outfile) } 
 
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
id2gofile=file.path(path.expand("~"),"data","mammalsid2go.tab") 
uniprot$id2go(file,id2gofile) 
head(read.table(id2gofile,sep="\t",quote="")) 
} 
 
