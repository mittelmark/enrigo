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
  the GO ids and the gene ids for a given species. For getting a list of available species and for the  
  download a internet connection is required. 
} 
\value{local filename, invisible} 
\examples{ 
# get all species 
fname=enr$download("Apis mellifera",folder=file.path(path.expand("~"),"data")) 
enr$new(fname) 
fname 
} 
 
