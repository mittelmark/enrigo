\name{enr$topplot} 
\alias{enr_topplot} 
\alias{enr$topplot} 
\title{Create an color coded plot of expression levels} 
\description{ 
  Creates a plot were RNA or other expression levels are encoded using more  
  or less saturated color schemas. 
} 
\usage{ enr_topplot(x, range=c(0,10),legend=TRUE,scale="red",text.start=-0.8) } 
 
\arguments{ 
 \item{x}{a two to six column data frame with expression levels, usually on a log2-scale. } 
 \item{range}{value range used for the scaling of colors} 
 \item{legend}{Should a color legend be drawn right of the color code rectangles, default: TRUE} 
 \item{scale}{The color schema used for scaling, possible values are 'red', 'green', 'darkgreen', 'blue' and 'gray' (or 'grey'), default: 'red'}  
 \item{text.start}{Where to place te rownames on the x-axis, change this if you have short labels to more positive values, default: -0.8} 
} 
\examples{ 
df=data.frame(matrix(runif(104,min=0,max=10),ncol=4)) 
colnames(df)=c("X1","X2", "X3", "X4") 
rownames(df)=LETTERS[1:26] 
df=df[order(df[,1]/df[,2]),] 
par(mfrow=c(1,2),mai=rep(0.1,4)) 
enr$topplot(df,text.start=0.8,scale="darkgreen") 
enr$topplot(df,text.start=0.8,scale="red") 
} 
 
