\name{goutil$download} 
\alias{goutil_download} 
\alias{goutil$download} 
\title{ Download for the given year the obofile } 
\description{ 
  This function is used to download the Obo file for the given year from the gene ontology side. 
} 
\usage{ goutil_download(version,timeout=600,folder=NULL,...) } 
 
\arguments{ 
  \item{version}{either number such as 2005 or higher or character string for a version name like "2023-01-01" } 
  \item{timeout}{allowed timeout in seconds how long the download might need, large downloads might take a long, time, default: 600} 
  \item{folder}{data folder where theGO obo data files should be stored, if not given, the current working directory is used, default: NULL} 
  \item{\ldots}{Other arguments delegated to the function download.file} 
} 
\details{ 
  This function downloads the GO obofile from the Gene Ontology website. 
} 
\value{filename of the downloaded file} 
\examples{ 
obofile=goutil$download(2022, 
folder=file.path(path.expand("~"),"data")) 
} 
 
