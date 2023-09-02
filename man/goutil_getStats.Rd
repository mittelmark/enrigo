\name{goutil$getStats} 
\alias{goutil_getStats} 
\alias{goutil$getStats} 
\title{Return statistics for the current GO file} 
\description{ 
This function returns a data frame with a overal summary statistics for the currently loaded GO file. 
} 
\usage{ goutil_getStats() } 
 
\value{data frame with the columns active and obsolete and the rownames c, f, p for the three namespaces} 
\examples{ 
goutil$getStats() 
} 
 
