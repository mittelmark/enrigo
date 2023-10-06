\name{enr$annotation} 
\alias{enr_annotation} 
\alias{enr$annotation} 
\title{Extract GO annotations for the given genes} 
\description{ 
  This function takes a vector of gene indentifiers and extracts for them the annotations 
  from the current GAF file. 
} 
\usage{ enr_annotation(x) } 
 
\arguments{ 
  \item{x}{vector of gene identifiers} 
} 
\details{ 
  This function allows you investigate for a given set of genes what are the available annotations. 
} 
\value{data frame with the columns from the annotation file, you can add yourself the GO names as shown in the example} 
\examples{ 
enr$gaf("Apis mellifera") 
anno=enr$annotation(c("LOC102653920","LOC725343","DUMMY")) 
head(anno) 
fname=goutil$download(2022) 
 
goutil$new(fname) 
anno=cbind(anno,go_name=goutil$getName(anno$go_id)) 
head(anno) 
} 
 
