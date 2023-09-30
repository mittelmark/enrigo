\name{uniprot$download} 
\alias{uniprot_download} 
\alias{uniprot$download} 
\title{ Download SwissProt UniProt data } 
\description{ 
  This function downloads SwissProt database files of the UniProt database 
  to extract from the file GO annotation data. 
} 
\usage{ uniprot_download(taxa,timeout=600,folder=NULL) } 
 
\arguments{ 
  \item{taxa}{ The taxonomic division, valid values are 
 "archaea","bacteria","fungi","human","invertebrates","mammals","rodents","vertebrates" or "viruses"} 
  \item{timeout}{ The allowed timeout in seconds how long the download might need, large downloads might take a long, time, default: 600} 
  \item{folder}{the data folder where the uniprot data files should be stored, if not given, the current working directory is used, default: NULL} 
} 
\details{ 
  This function allows you to prepare annotation file on your own if you do not can ulize the existing 
  annotation files from the GO consortium. To download the data 
  a internet connection is required. 
} 
\value{local filename, invisible} 
\examples{ 
# get all species 
file=uniprot$download("mammals", 
  folder=file.path(path.expand("~"),"data")) 
file 
} 
 
