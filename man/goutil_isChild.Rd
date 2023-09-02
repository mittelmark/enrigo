\name{goutil$isChild} 
\alias{goutil_isChild} 
\alias{goutil$isChild} 
\title{Checks if the given child is indeed a child node of parent} 
\description{ 
This function checks for a given child node if it is a  
child node of the parent node. 
} 
\usage{ goutil_isChild(parent,child) } 
 
\arguments{ 
  \item{parent}{tested parent GO id} 
  \item{child}{tested child GO id} 
} 
\value{boolean, TRUE if parent id is a parent node of the child id} 
\examples{ 
  goutil$new("2023-01-01-go.obo") 
  goutil$isChild("GO:0003674","GO:0005488") # should be TRUE 
  goutil$isChild("GO:0003674","GO:0000001") # should be FALSE 
} 
 
