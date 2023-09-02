\name{goutil$getOboVersion} 
\alias{goutil_getOboVersion} 
\alias{goutil$getOboVersion} 
\title{Returns the actual version date for the current obo file} 
\description{ 
This function returns a the editing date for the current OBO file in the form YYYY-MM-DD as string. 
} 
\usage{ goutil_getOboVersion() } 
 
\value{editing date of the Obo file as a string} 
\examples{ 
obofile=goutil$new("2023-01-01-go.obo") 
goutil$getOboVersion() 
} 
 
