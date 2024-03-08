\name{enr$enrichment} 
\alias{enr_enrichment} 
\alias{enr$enrichment} 
\title{Perform a GO term enrichment analyzes for the given gene lists} 
\description{ 
  This function takes two gene lists and performs a enrichment analysis. 
   
} 
\usage{ enr_enrichment(fullset,subset,mapping=NULL,max.genes=5, derichment=FALSE, min.term=2, namespace="all", alternative="greater")} 
 
\arguments{ 
  \item{fullset}{full set of genes} 
  \item{subset}{subset of genes, for instance these ones which are upregulated} 
  \item{mapping}{mapping file for any enrichment analysis, must contain a column id and a column gene, if not given GO enrichment is done} 
  \item{max.genes}{maximal number of gene identifiers attached into the result table, default: 5} 
  \item{derichment}{should as well a de-richment analysis performed, default: FALSE} 
  \item{min.term}{minimal number which should a GO id or other terms should be present, default: 2} 
  \item{namespace}{which namespace terms should be checked, 'c','f','p' or 'all', default: 'all'} 
  \item{alternative}{should the fisher.test performed one-sided only ("greater", lower p-values in some cases) as only enrichement is checked ) or 'two.sided', , default: 'greater'} 
} 
\details{ 
  This function performs the actual enrichment analysis. 
  In case the option `derichment` is set to TRUE, all genes of the fullset 
  will be checked for enrichment and derichment, if this is FALSE only the GO's 
  of the subset will be checked. We recommend the adjustment method 'BY' as it takes into consideration that the GO terms are not independent from each other. 
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
df=cbind(df,p_slim=tfslim) 
head(df,20) 
## term enrichment for your own terms 
set.seed(123) 
terms=sprintf("TRM\%03i",1:10) 
head(terms) 
genes=sprintf("GEN\%03i",1:50) 
head(genes) 
df=data.frame(gene=c(),term=c()) 
for (g in genes) {  
   x=sample(1:5,1);  
   df=rbind(df,data.frame(gene=rep(g,x),term=sample(terms,x))) 
 } 
head(df) 
# enrich term TRM001 for demo purposes 
idx=which(df[,'term'] == "TRM001") 
idx=idx[1:(length(idx)-2)]  
enr$enrichment(unique(df$gene), unique(df$gene[c(1:6,idx)]),df) 
} 
 
