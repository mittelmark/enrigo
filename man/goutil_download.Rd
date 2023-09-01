\name{goutil$download} 
\alias{goutil_download} 
\alias{goutil$download} 
\title{ Download for the given year the obofile } 
\description{ 
  This function is used to download the Obo file for the given year from the gene ontology side. 
} 
\usage{ goutil_download(version) } 
 
\arguments{ 
  \item{version}{ Either such as 2005 or higher or character string for a version name like "2023-01-01" } 
} 
\details{ 
  This function downloads the GO obofile from the Gene Ontolofy website. 
} 
\value{filename of the downloaded file} 
\examples{ 
obofile=goutil$download(2022) 
} 
 
