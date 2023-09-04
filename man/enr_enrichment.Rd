\name{enr$enrichment} 
\alias{enr_enrichment} 
\alias{enr$enrichment} 
\title{Perform a GO term enrichment analyzes for the given gene lists} 
\description{ 
  This function takes two gene lists and performs a enrichment analysis. 
   
} 
\usage{ enr_enrichment(fullset,subset) } 
 
\arguments{ 
  \item{fullset}{full set of genes} 
  \item{subset}{subset of genes, for instance these ones which are upregulated} 
} 
\details{ 
  This function performs the actual enrichment analysus. 
} 
\value{data frame with the following columns  
\itemize{ 
\item go_id - the GO id 
\item total - number of all genes from the given fullset found in the annotation file belonging to this GO id 
\item ntotal - number of all genes from the given fullset found in the annotation file annotation file not belonging to this GO id 
\item set - number of all genes from the given subset found in the annotation file belonging to this GO id 
\item nset - number of all genes from the given subset found in the annotation file not belonging to this GO id 
\item p_value - the raw p-value of the fisher test 
\item residuals - the pearson residuals for the contigency table 
\item cohens_w - the effect size Cohen-s w for the contingency table 
} 
} 
\examples{ 
fname=goutil$download(2022) 
goutil$new(fname) 
enr$gaf("Apis mellifera") 
enr$new() 
fullset=sample(unique(enr$data$godata$ID1),80) 
## some symbols 
fullset=c(fullset,sample(unique(enr$data$godata$ID2),10)) 
# add some non-annotated stupid gene ids 
fullset=c(fullset, sprintf("dummy\%02i",1:10)) 
# just select some random genes for enrichment analysis 
subset=sample(fullset,30)  
# How many genes of the subset are in the GAF file 
gaf.subset=intersect(enr$symbol2loc(subset),enr$data$godata$ID1) 
length(gaf.subset) 
# How many genes of the fullset are in the GAF file 
gaf.fullset=intersect(enr$symbol2loc(fullset),enr$data$godata$ID1) 
length(gaf.fullset) 
df=enr$enrichment(gaf.fullset,gaf.subset) 
head(df) 
dim(df) 
# Now let's only display the drosophila slims 
slims=goutil$getSlims(slim="drosophila",nsp="p") 
length(slims) 
tfslim=rep(FALSE,nrow(df)) 
tfslim[df$go_id \%in\% slims]=TRUE 
table(tfslim) 
enr=cbind(df,p_slim=tfslim) 
head(df,20) 
} 
 
