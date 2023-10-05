\name{enr$name2id} 
\alias{enr_name2id} 
\alias{enr$name2id} 
\title{Convert names into identifiers} 
\description{ 
  This function takes a vector of text strings, may be as well with duplicates 
  and converts them into identifiers.#'   identifiers and returns a vector of gene indentifiers based on the GO 
} 
\usage{ enr_name2id(x,pattern="TM\%04i") } 
 
\arguments{ 
  \item{x}{vector of text strings} 
  \item{pattern}{the id pattern which should be used to create the identifiers using sprintf, default: 'TM\%04i'} 
} 
\details{ 
  This function allows you to create unique identifiers for text strings. 
} 
\value{data frame with the columns id and and name, where the latter are the given text strings} 
\examples{ 
# get all species 
enr$name2id(c("hello","world","dummy","hello")) 
} 
 
