\name{enr$gaf} 
\alias{enr_gaf} 
\alias{enr$gaf} 
\title{ Setup or query a Gene Ontology annotation file } 
\description{ 
  This function setups the Gene Ontology Annotation file.  
   
} 
\usage{ enr_gaf(x=NULL,folder=NULL) } 
 
\arguments{ 
  \item{x}{ The filename or a species name for the file which should be used to perform the  
   enrichment analysis, if the file does not exists it can be download, if no filename is given all available species are listed, default: NULL } 
  \item{folder}{the data folder where the annotation files should be stored, if not given, the current working directory is used, default: NULL} 
} 
\details{ 
  This function allows you to intialize your annotation file which contains the mapping between 
  the GO ids and the gene ids for a given species. For getting a list of available species and for the  
  download a internet connection is required. 
} 
\value{local filename, invisible} 
\examples{ 
# get all species 
res=enr$gaf() 
sort(names(res)) 
enr$gaf("Apis mellifera",folder=file.path(path.expand("~"),"data")) 
} 
 
