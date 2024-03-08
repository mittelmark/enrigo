\name{enr$vennplot} 
\alias{enr_vennplot} 
\alias{enr$vennplot} 
\title{Create a two-set venn diagram with the standard color codes} 
\description{ 
  Use a TRUE/FALSE data frame to create a two set venn diagram plot. 
} 
\usage{ enr_vennplot(x, cols=c("#EF536BBB","#61E04FBB","#536BEFBB")) } 
 
\arguments{ 
 \item{x}{a two column table with TRUE/FALSE values} 
 \item{cols}{list of two or three colors with transparency,  
             default: c("#EF536BBB","#61E04FBB","#536BEFBB")} 
} 
\examples{ 
tab=data.frame(matrix(rnorm(100),ncol=2)) 
colnames(tab)=c('A', 'B') 
tab=tab>0 
head(tab) 
table(tab[,1],tab[,2]) 
enr$vennplot(tab) 
} 
 
