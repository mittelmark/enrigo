\name{enr$enrichment_ensembl} 
\alias{enr_enrichment_ensembl} 
\alias{enr$enrichment_ensembl} 
\title{Perform a enrichment analyzes for the given subset gene lists based on Ensembl genes} 
\description{ 
  This function takes a list of Ensembl gene identifiers which are checked for possible term enrichments. 
  The fullset is here taken from the Uniprot/Swissprot database 
} 
\usage{ enr_enrichment_ensembl(subset, max.genes=5, mapping=NULL, folder=file.path(path.expand("~"),"data")) } 
 
\arguments{ 
  \item{subset}{subset of genes, for instance these ones which are upregulated in case of non-human 
   genes you must submit human Ensembl genes homologous to them} 
  \item{max.genes}{maximal number of gene identifiers attached into the result table, default: 5} 
  \item{mapping}{mapping file for any enrichment analysis, must contain a column id and a column gene, if not given GO enrichment is done} 
  \item{folder}{the data folder where to store intermediate results} 
} 
\details{ 
  This function performs an enrichment analysis based on human Ensembl gene identifiers, 
  if you have non-human genes you should submit the human genes which are homologous to your species genes. 
} 
\value{data frame with the following columns  
\itemize{ 
\item go_id - the GO or term id  
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
genes=read.table(text="ENSG00000102057 
ENSG00000258659 
ENSG00000184206 
ENSG00000184345 
ENSG00000153684 
ENSG00000278662 
ENSG00000183629 
ENSG00000125520") 
head(genes[,1]) 
enr$enrichment_ensembl(genes[,1]) 
} 
 
