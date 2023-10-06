\name{enr$download} 
\alias{enr_download} 
\alias{enr$download} 
\title{Download a Gene Ontology annotation file for the given species} 
\description{ 
  This function downloads a Gene Ontology Annotation file.  
} 
\usage{ enr_download(x,folder=NULL) } 
 
\arguments{ 
  \item{x}{species name for the file which should be used to perform the  
   enrichment analysis, if the file does not already on your harddisk or if the file is not from the current month it is be download, default: NULL 
  if no filename is given all available species are listed, default: NULL } 
  \item{folder}{the data folder where the annotation files should be stored, if not given, the current working directory is used, default: NULL} 
} 
\details{ 
  This function allows you to download your annotation file which contains the mapping between 
  the GO ids and the gene ids for a given species. 
  This function together with the function enr$new allows you to do the same thing as with the enr$gaf  
  function but the approach first, download and get a filename then intialize with this file is probably more transparent for the user. 
} 
\value{local filename} 
\examples{ 
# get all species 
fname=enr$download("Apis mellifera",folder=file.path(path.expand("~"),"data")) 
fname 
enr$new(fname) 
} 
 
