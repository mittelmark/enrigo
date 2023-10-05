#' \docType{class}  
#' \name{enr} 
#' \alias{enr}
#' \alias{enr-class}
#' \title{ Environment object with functions do do enrichment analysis for GO terms }
#' \description{
#' The functions within the enr environment allow to do Gene Ontology term enrichment analysis for a given list
#' of genes against a full list of genes from the same organism.
#' }
#' \section{Methods}{
#' \itemize{
#' \item \code{\link[enrigo:enr_gaf]{enr$gaf(filename)}}{query, download or initialize a Gene Ontology association file}
#' }
#' }
#' \examples{
#' # get all species
#' res=enr$gaf()
#' sort(names(res))
#' } 

enr=new.env()


#' \name{enr$gaf}
#' \alias{enr_gaf}
#' \alias{enr$gaf}
#' \title{ Setup or query a Gene Ontology annotation file }
#' \description{
#'   This function setups the Gene Ontology Annotation file. 
#'   
#' }
#' \usage{ enr_gaf(x=NULL) }
#'
#' \arguments{
#'   \item{x}{ The filename or a species name for the file which should be used to perform the 
#'    enrichment analysis, if the file does not exists it can be download, if no filename is given all available species are listed, default: NULL }
#' }
#' \details{
#'   This function allows you to intialize your annotation file which contains the mapping between
#'   the GO ids and the gene ids for a given species. For getting a list of available species and for the 
#'   download a internet connection is required.
#' }
#' \value{local filename, invisible}
#' \examples{
#' # get all species
#' res=enr$gaf()
#' sort(names(res))
#' }
#' 

enr$gaf = function (x=NULL) {
    if (is.null(x)) {
        # check which species are available
        download.file("http://current.geneontology.org/products/pages/downloads.html","download.html")
        fin  = file("download.html", "r")
        spec=list()
        while(length((line = readLines(fin,n=1)))>0) {
            if (grepl("<b>[A-Z][a-z]+ [a-z]+</b>",line)) {
                s=gsub(".+<b>([A-Z][a-z]+ [a-z]+)</b>","\\1",line)
            } else if (grepl("<a href=\"http://current.geneontology.org/annotations/.+?gaf.gz",line)) {
                f=gsub(".+annotations/(.+?.gaf.gz).+","\\1",line)
                f=paste("http://current.geneontology.org/annotations/",f,sep="")
                spec[[s]]=f
            }
        }
        close(fin)
        download.file("https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/","download.html")
        fin  = file("download.html", "r")
        while(length((line = readLines(fin,n=1)))>0) {
            if (grepl("<a href=\".+HGD_go_annotation.gaf.gz",line)) {
                s=gsub("_"," ",gsub(".+<a href=\"(.+?)_HGD.+","\\1",line))
                f=gsub(".+<a href=\"(.+?_HGD_go_annotation.gaf.gz)\">.+","\\1",line)
                f=paste("https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/",f,sep="")
                spec[[s]]=f
            }
        }
        close(fin)
        return(spec)
    } else if (!grepl(".gaf",x) & !grepl(".gaf.gz",x)) {
        # must be a species
        spec=enr$gaf()
        if (!x %in% names(spec)) {
            stop("Error: Given species is not available!\nValid species names are:", paste(sort(names(spec)),sep="',\n'"))
        }
        dfile=gsub(" ","_",paste(x,".gaf.gz",sep=""))
        download.file(spec[[x]], destfile=dfile)
        enr$annotation=dfile
    } else  if (!file.exists(x)) {
        stop(paste("Error: file",x,"does not exists!"))
    } else if (grepl(".gaf",x)) {
        enr$annotation=x
        invisible(enr$annotion)
    } else {
        stop("Invalid filename or species, filenames must end with gaf or gaf.gz")
    }
}

#' \name{enr$new}
#' \alias{enr_new}
#' \alias{enr$new}
#' \title{Initialize the annotation data }
#' \description{
#'   This function reads in the Gene Ontology annotation file. 
#'   
#' }
#' \usage{ enr_new(gaffile=NULL) }
#'
#' \arguments{
#'   \item{gaffile}{GO annotation filename which should be used to perform the 
#'    enrichment analysis, if not given the last downloaded file, using the function `enr$gaf` is used, default: NULL }
#' }
#' \details{
#'   This function allows you to initialize the annotation data which contain the mapping between
#'   the GO ids and the gene ids for a given species.
#' }
#' \value{NULL}
#' \examples{
#' # get all species
#' enr$gaf("Apis mellifera")
#' enr$new()
#' }
#' 

enr$new <- function (gaffile=NULL) {
    self=enr
    if (is.null(gaffile)) {
        gaffile=self$annotation
    }
    godata=read.table(gaffile,comment.char="!",sep="\t",quote="")
    godata=godata[,c(2,3,5,7)]
    colnames(godata)=c("ID1","ID2","GOID","EvCode")
    if (grepl("Apis.+mellif",gaffile)) {
        godata$ID1=paste("LOC",godata$ID1,sep="")
    }
    self$data=list()
    self$data$godata=godata
}

#' \name{enr$enrichment}
#' \alias{enr_enrichment}
#' \alias{enr$enrichment}
#' \title{Perform a GO term enrichment analyzes for the given gene lists}
#' \description{
#'   This function takes two gene lists and performs a enrichment analysis.
#'   
#' }
#' \usage{ enr_enrichment(fullset,subset,mapping=NULL) }
#'
#' \arguments{
#'   \item{fullset}{full set of genes}
#'   \item{subset}{subset of genes, for instance these ones which are upregulated}
#'   \item{mapping}{mapping file for any enrichment analysis, must contain a column id and a column gene, if not given GO enrichment is done}
#' }
#' \details{
#'   This function performs the actual enrichment analysus.
#' }
#' \value{data frame with the following columns 
#' \itemize{
#' \item go_id - the GO id 
#' \item total - number of all genes from the given fullset found in the annotation file belonging to this GO id
#' \item ntotal - number of all genes from the given fullset found in the annotation file annotation file not belonging to this GO id
#' \item set - number of all genes from the given subset found in the annotation file belonging to this GO id
#' \item nset - number of all genes from the given subset found in the annotation file not belonging to this GO id
#' \item p_value - the raw p-value of the fisher test
#' \item residuals - the pearson residuals for the contigency table
#' \item cohens_w - the effect size Cohen-s w for the contingency table
#' }
#' }
#' \examples{
#' fname=goutil$download(2022)
#' goutil$new(fname)
#' enr$gaf("Apis mellifera")
#' enr$new()
#' fullset=sample(unique(enr$data$godata$ID1),80)
#' ## some symbols
#' fullset=c(fullset,sample(unique(enr$data$godata$ID2),10))
#' # add some non-annotated stupid gene ids
#' fullset=c(fullset, sprintf("dummy\%02i",1:10))
#' # just select some random genes for enrichment analysis
#' subset=sample(fullset,30) 
#' # How many genes of the subset are in the GAF file
#' gaf.subset=intersect(enr$symbol2loc(subset),enr$data$godata$ID1)
#' length(gaf.subset)
#' # How many genes of the fullset are in the GAF file
#' gaf.fullset=intersect(enr$symbol2loc(fullset),enr$data$godata$ID1)
#' length(gaf.fullset)
#' df=enr$enrichment(gaf.fullset,gaf.subset)
#' head(df)
#' dim(df)
#' # Now let's only display the drosophila slims
#' slims=goutil$getSlims(slim="drosophila",nsp="p")
#' length(slims)
#' tfslim=rep(FALSE,nrow(df))
#' tfslim[df$go_id \%in\% slims]=TRUE
#' table(tfslim)
#' enr=cbind(df,p_slim=tfslim)
#' head(df,20)
#' ## term enrichment for your own terms
#' set.seed(123)
#' terms=sprintf("TRM\%03i",1:10)
#' head(terms)
#' genes=sprintf("GEN\%03i",1:50)
#' head(genes)
#' df=data.frame(gene=c(),term=c())
#' for (g in genes) { 
#'    x=sample(1:5,1); 
#'    df=rbind(df,data.frame(gene=rep(g,x),term=sample(terms,x)))
#'  }
#' head(df)
#' # enrich term TRM001 for demo purposes
#' idx=which(df[,'term'] == "TRM001")
#' idx=idx[1:(length(idx)-2)] 
#' ##enr$enrichment(unique(df$gene), unique(df$gene[c(1:6,idx)]),df)
#' }
#' 

enr$enrichment <- function (fullset,subset,mapping=NULL) {
    cohensW = function (x,p=NULL) {
        if (is.table(x) | is.matrix(x)) {
            tab=x
            pe=prop.table(chisq.test(tab)$expected)
            po=prop.table(tab)
            w=sqrt(sum(((po-pe)^2)/pe))
            return(w[[1]])
        } else if (is.null(p)) {
            stop('Error: If x is a vector, p must be given!')
        } else {
            if (length(x) == 2 & length(p) == 1) {
                p=c(p,1-p)
                po=prop.table(x)
                pe=p
            } else if  (length(x) == length(p)) {
                po=prop.table(x)
                pe=p
            } else {
                stop('Error: for more than 2 categories the
                     given proportion vector p must have the
                     same length as the given count vector x')
            }
            w = sqrt(sum(((po-pe)^2)/pe))
            return(w)
        }
    }
    df = data.frame(go_id = as.character(c()),
                    total  = as.numeric(c()),
                    ntotal  = as.numeric(c()),
                    set    = as.numeric(c()),
                    nset    = as.numeric(c()),
                    p_value = as.numeric(c()),
                    residuals = as.numeric(c()),
                    cohens_w = as.numeric(c())
                    )
    if (is.data.frame(mapping) | is.matrix(mapping)) {
        # total number of genes
        # some genes might be not in mapping
        fmapping = mapping[mapping[,'gene'] %in% fullset,]
        smapping = mapping[mapping[,'gene'] %in% subset,]
        CA = length(unique(fmapping[,'gene']))
        DA = length(unique(smapping[,'gene']))
        owarn=options("warn")[[1]]
        options(warn=-1)
        for (term in unique(fmapping[,'term'])) {
            # A how many genes in the fullset have this annotation
            A = length(unique(fmapping[fmapping[,'term'] == term,'gene']))
            # B how many genes in the subset have this annotation
            B = length(unique(smapping[smapping[,'term'] == term,'gene']))
            C=CA-A
            D=DA-B 
            #how many genes are not annotated with this GO in the gaf file
            mat=matrix(c(A,B,C,D) , ncol=2,byrow=TRUE)
            
            pval = fisher.test(mat)$p.value
            res = chisq.test(mat)$residuals[1,2]
            cohW = cohensW(mat)
            total = A
            set   = B
            df = rbind(df, 
                       data.frame(
                                  go_id = term,
                                  total  = A,
                                  ntotal = C,
                                  set    = B,
                                  nset   = D,
                                  p_value = pval,
                                  residuals = res,
                                  cohens_w = cohW)
                       )
        }
            
    } else {
        fullset=unique(enr$symbol2loc(fullset))
        idx1=which(enr$data$godata$ID1 %in% fullset)
        subset=unique(enr$symbol2loc(subset))
        godata1=enr$data$godata[idx1,c('ID1','GOID')]
        idx3=which(enr$data$godata$ID1 %in% subset)
        #print(length(idx))
        godata2=enr$data$godata[idx3,c('ID1','GOID')]
        godata1$GOID = goutil$altid2new(godata1$GOID)
        godata2$GOID = goutil$altid2new(godata2$GOID)
        CA = length(unique(godata1$ID1))
        #how many go ids are in the list
        DA = length(unique(godata2$ID1))
        owarn=options("warn")[[1]]
        options(warn=-1)
        for (go in unique(godata2$GOID)) {
            #how many genes have this GO annotation is in the file
            A = length(unique(godata1[godata1$GOID == go,'ID1']))
            #how many times, this GO is in the given geneIds list
            B = length(unique(godata2[godata2$GOID == go,'ID1']))
            C=CA-A
            D=DA-B 
            #how many genes are not annotated with this GO in the gaf file
            mat=matrix(c(A,B,C,D) , ncol=2,byrow=TRUE)
            
            pval = fisher.test(mat)$p.value
            res = chisq.test(mat)$residuals[1,2]
            cohW = cohensW(mat)
            total = A
            set   = B
            df = rbind(df, 
                       data.frame(
                                  go_id = go,
                                  total  = A,
                                  ntotal = C,
                                  set    = B,
                                  nset   = D,
                                  p_value = pval,
                                  residuals = res,
                                  cohens_w = cohW)
                       )
        }
        names     = goutil$getName(df$go_id)
        df=cbind(df,go_name=names)
    }
    # sort the data frame base on the cohens_w
    df = df[order(df$cohens_w, decreasing = TRUE),]
    df=cbind(df,fdr=p.adjust(df$p_value,method="BH"))
    options(warn=owarn)
    return(df)
}  

#' \name{enr$symbol2loc}
#' \alias{enr_symbol2loc}
#' \alias{enr$symbol2loc}
#' \title{Convert possible gene symbols into database identifiers}
#' \description{
#'   This function takes a vector of gene symbols or gene symbols mixed with gene
#'   identifiers and returns a vector of gene indentifiers based on the GO
#'   annotation file.
#' }
#' \usage{ enr_symbol2loc(symbol) }
#'
#' \arguments{
#'   \item{symbol}{vector of gene symbols or gene symbols mixed with gene identifiers}
#' }
#' \details{
#'   This function allows you select the mapping gene identifiers if you have only gene symbols in your gene list.
#' }
#' \value{vector of NCBI gene identifiers, LOC values}
#' \examples{
#' # get all species
#' enr$gaf("Apis mellifera")
#' enr$new()
#' enr$symbol2loc(c('Lys-3','Dbp80',
#'   'LOC727025','LOC102653920','Dummy1'))
#' }
#' 

enr$symbol2loc <- function (symbol) {
    self=enr
    if (length(symbol)>1) {
        return(unlist(lapply(symbol,function (x) { return(enr$symbol2loc(x)) })))
    } else {
        idx = which(enr$data$godata$ID2 == symbol)
        if (length(idx)>0) {
            return(enr$data$godata$ID1[idx[1]])
        } else {
            return(symbol)
        }
    }
}

#' \name{enr$name2id}
#' \alias{enr_name2id}
#' \alias{enr$name2id}
#' \title{Convert names into identifiers}
#' \description{
#'   This function takes a vector of text strings, may be as well with duplicates
#'   and converts them into identifiers.#'   identifiers and returns a vector of gene indentifiers based on the GO
#' }
#' \usage{ enr_name2id(x,pattern="TM\%04i") }
#'
#' \arguments{
#'   \item{x}{vector of text strings}
#'   \item{pattern}{the id pattern which should be used to create the identifiers using sprintf, default: 'TM\%04i'}
#' }
#' \details{
#'   This function allows you to create unique identifiers for text strings.
#' }
#' \value{data frame with the columns id and and name, where the latter are the given text strings}
#' \examples{
#' # get all species
#' enr$name2id(c("hello","world","dummy","hello"))
#' }
#' 

enr$name2id <- function (x,pattern="TM%04i") {
    k=1
    res=list()
    df=data.frame(term=c(),name=c())
    for (i in x) {
        if (is.null(res[[i]])) {
            term=sprintf(pattern,k)
            res[[i]]=term
            k=k+1
        } else {
            term=res[[i]]
        }
        df=rbind(df,data.frame(id=term,name=i))
    }
    return(df)
}
enr_enrichment = enr$enrichment
enr_gaf        = enr$gaf
enr_new        = enr$new
enr_symbol2loc = enr$symbol2loc
enr_name2id    = enr$name2id
