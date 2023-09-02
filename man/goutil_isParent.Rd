\name{goutil$isParent} 
\alias{goutil_isParent} 
\alias{goutil$isParent} 
\title{Checks if the given parent node is indeed a parent node of child} 
\description{ 
This function checks for a given parent node if it is a  
parent node of the child node. 
} 
\usage{ goutil_isParent(parent,child) } 
 
\arguments{ 
  \item{parent}{tested parent GO id} 
  \item{child}{tested child GO id} 
} 
\value{boolean, TRUE if parent id is a parent node of the child id} 
\examples{ 
  goutil$new("2023-01-01-go.obo") 
  goutil$isParent("GO:0003674","GO:0005488") # should be TRUE 
  goutil$isParent("GO:0003674","GO:0000001") # should be FALSE 
} 
 
