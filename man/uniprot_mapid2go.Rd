\name{uniprot$mapid2go} 
\alias{uniprot_mapid2go} 
\alias{uniprot$mapid2go} 
\title{Create a mapping between your identifiers and the GO identifiers using identifiers provided by UniProt files} 
\description{ 
  This function creates a mapping file between your identifiers using 
  Uniprot identifiers to the GO identifiers provided by the SwissProt database. 
} 
\usage{ uniprot_mapid2go(upid2gofile,idfile,outfile,evidence="ALL") } 
 
\arguments{ 
  \item{upid2gofile}{file created with `uniprot$id2go`} 
  \item{idfile}{mapping file for instance created with `uniprot$ensg2go`} 
  \item{outfile}{outfile name where the mapping is stored, a tabfile} 
  \item{evidence}{The evidence codes which should be used, currently not used, default: "ALL"} 
} 
\value{NULL} 
\examples{ 
# get all species 
file=uniprot$download("mammals", 
  folder=file.path(path.expand("~"),"data")) 
id2gofile=file.path(path.expand("~"),"data","mammalsid2go.tab") 
uniprot$id2go(file,id2gofile) 
id2ensgfile=file.path(path.expand("~"),"data","id2ensg.tab") 
uniprot$id2ensg(file,id2ensgfile) 
mappingfile=file.path(path.expand("~"),"data","mapping.tab") 
uniprot$mapid2go(id2gofile,id2ensgfile,mappingfile) 
head(read.table(mappingfile,sep="\t",quote="")) 
} 
 
