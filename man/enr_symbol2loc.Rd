\name{enr$symbol2loc} 
\alias{enr_symbol2loc} 
\alias{enr$symbol2loc} 
\title{Convert possible gene symbols into database identifiers} 
\description{ 
  This function takes a vector of gene symbols or gene symbols mixed with gene 
  identifiers and returns a vector of gene indentifiers based on the GO 
  annotation file. 
} 
\usage{ enr_symbol2loc(symbol) } 
 
\arguments{ 
  \item{symbol}{vector of gene symbols or gene symbols mixed with gene identifiers} 
} 
\details{ 
  This function allows you select the mapping gene identifiers if you have only gene symbols in your gene list. 
} 
\value{vector of NCBI gene identifiers, LOC values} 
\examples{ 
# get all species 
enr$gaf("Apis mellifera") 
enr$new() 
enr$symbol2loc(c('Lys-3','Dbp80', 
  'LOC727025','LOC102653920','Dummy1')) 
} 
 
