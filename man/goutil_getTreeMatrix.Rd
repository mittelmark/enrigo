\name{goutil$getTreeMatrix} 
\alias{goutil_getTreeMatrix} 
\alias{goutil$getTreeMatrix} 
\title{Return an adjacency matrix for the given ids in tree} 
\description{ 
 This function returns an adjacency matrix for the given list of 
 GO ids in the form of an adjacency matrix. 
} 
\usage{ goutil_getTreeMatrix(tree) } 
 
\arguments{ 
  \item{tree}{a list of GO ids, usually created with the goutil$getTee function or a single GO id, in this case the `getTree` function will be called internally} 
} 
\examples{ 
 goutil$new("2023-01-01-go.obo") 
 tree=goutil$getTree("GO:0003676") 
 goutil$getTreeMatrix(tree) 
} 
 
