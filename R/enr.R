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
#' \item \code{\link[enrigo:enr_annotation]{enr$annotation(genes)}}{return GO annotations for the given gene identifiers}
#' \item \code{\link[enrigo:enr_enrichment]{enr$enrichment(fullset,subset)}}{perform an enrichment analysis}
#' \item \code{\link[enrigo:enr_gaf]{enr$gaf()}}{query, download or initialize a Gene Ontology association file}
#' \item \code{\link[enrigo:enr_download]{enr$download(species)}}{dowmload an annotation file for the given species}
#' \item \code{\link[enrigo:enr_new]{enr$new(filename)}}{initialize the GO annotation}
#' \item \code{\link[enrigo:enr_name2id]{enr$name2id(names)}}{convert free text terms into identifiers}
#' \item \code{\link[enrigo:enr_symbol2loc]{enr$symbol2loc(symbol)}}{convert gene esymbols into LOC identifiers}
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
#' \usage{ enr_gaf(x=NULL,folder=NULL) }
#'
#' \arguments{
#'   \item{x}{ The filename or a species name for the file which should be used to perform the 
#'    enrichment analysis, if the file does not exists it can be download, if no filename is given all available species are listed, default: NULL }
#'   \item{folder}{the data folder where the annotation files should be stored, if not given, the current working directory is used, default: NULL}
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
#' enr$gaf("Apis mellifera",folder=file.path(path.expand("~"),"data"))
#' }
#' 

enr$gaf = function (x=NULL,folder=NULL) {
    if (is.null(x)) {
        # check which species are available
        if (!file.exists("download-go.html") || substr(file.mtime("download-go.html"),1,10) != substr(Sys.Date(),1,10)) {
            download.file("http://current.geneontology.org/products/pages/downloads.html","download-go.html")
        } 
        fin  = file("download-go.html", "r")
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
        if (!file.exists("download-hgd.html") || substr(file.mtime("download-hgd.html"),1,10) != substr(Sys.Date(),1,10)) {
            download.file("https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/","download-hgd.html")
        }
        fin  = file("download-hgd.html", "r")
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
        dfile=enr$download(x)
        enr$new()
        invisible(enr$.annotation)
    } else  if (!file.exists(x)) {
        stop(paste("Error: file",x,"does not exists!"))
    } else if (grepl(".gaf$",x) | grepl(".gaf.gz$",x)) {
        enr$.annotation=x
        enr$new()
        invisible(enr$.annotation)
    } else {
        stop("Invalid filename or species, filenames must end with gaf or gaf.gz")
    }
}

#' \name{enr$download}
#' \alias{enr_download}
#' \alias{enr$download}
#' \title{Download a Gene Ontology annotation file for the given species}
#' \description{
#'   This function downloads a Gene Ontology Annotation file. 
#' }
#' \usage{ enr_download(x,folder=NULL) }
#'
#' \arguments{
#'   \item{x}{species name for the file which should be used to perform the 
#'    enrichment analysis, if the file does not already on your harddisk or if the file is not from the current month it is be download, default: NULL
#'   if no filename is given all available species are listed, default: NULL }
#'   \item{folder}{the data folder where the annotation files should be stored, if not given, the current working directory is used, default: NULL}
#' }
#' \details{
#'   This function allows you to download your annotation file which contains the mapping between
#'   the GO ids and the gene ids for a given species.
#'   This function together with the function enr$new allows you to do the same thing as with the enr$gaf 
#'   function but the approach first, download and get a filename then intialize with this file is probably more transparent for the user.
#' }
#' \value{local filename}
#' \examples{
#' # get all species
#' fname=enr$download("Apis mellifera",folder=file.path(path.expand("~"),"data"))
#' fname
#' enr$new(fname)
#' }
#' 

enr$download <- function (x,folder=NULL) {
    spec=enr$gaf()
    if (!x %in% names(spec)) {
        stop("Error: Given species is not available!\nValid species names are:", paste(sort(names(spec)),sep="',\n'"))
    }
    dfile=gsub(" ","_",paste(x,".gaf.gz",sep=""))
    if (!is.null(folder)) {
        if (!dir.exists(folder)) {
            dir.create(folder,recursive=TRUE)
        }
        dfile=file.path(folder,dfile)
    } else {
        dfile=file.path(getwd(),dfile)
    } 
    if (!file.exists(dfile) || substr(file.mtime(dfile),1,7) != substr(Sys.Date(),1,7)) {
        download.file(spec[[x]], destfile=dfile)
    }
    enr$.annotation=dfile
    return(dfile)
}

#' \name{enr$new}
#' \alias{enr_new}
#' \alias{enr$new}
#' \title{Initialize the annotation data for a given file}
#' \description{
#'   This function reads in a Gene Ontology annotation file. 
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
#' apisfile=enr$gaf("Apis mellifera",folder=file.path(path.expand("~"),"data"))
#' cerafile=enr$gaf("Apis cerana",folder=file.path(path.expand("~"),"data"))
#' enr$new(cerafile)
#' }
#' 

enr$new <- function (gaffile=NULL) {
    self=enr
    if (is.null(gaffile)) {
        gaffile=self$.annotation
    }
    godata=read.table(gaffile,comment.char="!",sep="\t",quote="")
    godata=godata[,c(2,3,5,7)]
    colnames(godata)=c("ID1","ID2","GOID","EvCode")
    ## TODO: check for LOC identifiers and add them if neccessary
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
#' \usage{ enr_enrichment(fullset,subset,mapping=NULL,max.genes=5, derichment=FALSE, min.term=2, namespace="all", alternative="greater")}
#'
#' \arguments{
#'   \item{fullset}{full set of genes}
#'   \item{subset}{subset of genes, for instance these ones which are upregulated}
#'   \item{mapping}{mapping file for any enrichment analysis, must contain a column id and a column gene, if not given GO enrichment is done}
#'   \item{max.genes}{maximal number of gene identifiers attached into the result table, default: 5}
#'   \item{derichment}{should as well a de-richment analysis performed, default: FALSE}
#'   \item{min.term}{minimal number which should a GO id or other terms should be present, default: 2}
#'   \item{namespace}{which namespace terms should be checked, 'c','f','p' or 'all', default: 'all'}
#'   \item{alternative}{should the fisher.test performed one-sided only ("greater", lower p-values in some cases) as only enrichement is checked ) or 'two.sided', , default: 'greater'}
#' }
#' \details{
#'   This function performs the actual enrichment analysis.
#'   In case the option `derichment` is set to TRUE, all genes of the fullset
#'   will be checked for enrichment and derichment, if this is FALSE only the GO's
#'   of the subset will be checked. We recommend the adjustment method 'BY' as it takes into consideration that the GO terms are not independent from each other.
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
#' df=cbind(df,p_slim=tfslim)
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
#' enr$enrichment(unique(df$gene), unique(df$gene[c(1:6,idx)]),df)
#' }
#' 

enr$enrichment <- function (fullset,subset,mapping=NULL,max.genes=5,derichment=FALSE,min.term=2,namespace="all",alternative='greater') {
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
                    cohens_w = as.numeric(c()),
                    genes    = as.character(c()))
    alt=alternative

    if (is.data.frame(mapping) | is.matrix(mapping)) {
        # total number of genes
        # some genes might be not in mapping
        colnames(mapping)=gsub("_id","",colnames(mapping))
        fmapping = mapping[mapping[,'gene'] %in% fullset,]
        smapping = mapping[mapping[,'gene'] %in% subset,]
        CA = length(unique(fmapping[,'gene']))
        DA = length(unique(smapping[,'gene']))
        owarn=options("warn")[[1]]
        options(warn=-1)
        df=list()
        if (derichment) {
            terms = unique(fmapping[,'term'])
        } else {
            terms = unique(smapping[,'term'])
        }
        sterms=names(table(terms)>min.term)
        terms=sterms
        for (term in terms) {
            # A how many genes in the fullset have this annotation
            A = length(unique(fmapping[fmapping[,'term'] == term,'gene']))
            # B how many genes in the subset have this annotation
            genes=unique(smapping[smapping[,'term'] == term,'gene'])
            B = length(genes)
            C=CA-A
            D=DA-B 
            #how many genes are not annotated with this GO in the gaf file
            mat=matrix(c(C,A,D,B) , ncol=2,byrow=TRUE)
            #,
            pval = fisher.test(mat,alternative=alt)$p.value
            res = chisq.test(mat)$residuals[2,2]
            cohW = cohensW(mat)
            total = A
            set   = B
            if (B > max.genes) {
                genes=paste(paste(genes[1:max.genes],collapse=", "),", ...")
            } else {
                genes=paste(paste(genes,collapse=", "))
            }
            df[[term]]=data.frame(
                                  go_id = term,
                                  total  = A,
                                  ntotal = C,
                                  set    = B,
                                  nset   = D,
                                  p_value = pval,
                                  residuals = res,
                                  cohens_w = cohW,
                                  genes=genes)
         }
         df=do.call(rbind,df)
            
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
        df=list()
        if (derichment) {
            terms = unique(godata1$GOID)
        } else {
            terms = unique(godata2$GOID)
        }
        if(namespace %in% c('c','p','f')) {
            nsp=goutil$getNamespace(terms)
            terms=terms[nsp == namespace]
        } else if (!namespace == "all") {
            stop("Error: Wrong namespace,valid names are 'c','p','f' or 'all'!")
        }
        for (go in terms) {
            #how many genes have this GO annotation is in the file
            A = length(unique(godata1[godata1$GOID == go,'ID1']))
            #how many times, this GO is in the given geneIds list
            genes=unique(godata2[godata2$GOID == go,'ID1'])
            B = length(genes)
            if (B > max.genes) {
                genes=paste(paste(genes[1:max.genes],collapse=", "),", ...")
            } else {
                genes=paste(paste(genes,collapse=", "))
            }

            C=CA-A
            D=DA-B 
            #how many genes are not annotated with this GO in the gaf file
            mat=matrix(c(C,A,D,B) , ncol=2,byrow=TRUE)
            pval = fisher.test(mat,alternative=alt)$p.value
            res = chisq.test(mat)$residuals[2,2]
            cohW = cohensW(mat)
            total = A
            set   = B
            df[[go]] = data.frame(
                                  go_id = go,
                                  total  = A,
                                  ntotal = C,
                                  set    = B,
                                  nset   = D,
                                  p_value = pval,
                                  residuals = res,
                                  cohens_w = cohW,
                                  genes=genes)

            #df = rbind(df, )
        }
        df=do.call(rbind,df)
        names     = goutil$getName(df$go_id)
        df=cbind(df,go_name=names)
        df=cbind(df,go_nsp=goutil$getNamespace(df$go_id))
    }
    # sort the data frame base on the cohens_w
    df = df[order(df$cohens_w, decreasing = TRUE),]
    if (min.term>1) {
        df=df[df$set>=min.term,]
    } 
    if (namespace !="all") {
        df=df[df$go_nsp==namespace,]
    }
    # exlude negative residuals
    df=df[df$residuals>0,]
    df=cbind(df,fdr=p.adjust(df$p_value,method='BH'))
    options(warn=owarn)
    return(df)
}  

#' \name{enr$enrichment_ensembl}
#' \alias{enr_enrichment_ensembl}
#' \alias{enr$enrichment_ensembl}
#' \title{Perform a enrichment analyzes for the given subset gene lists based on Ensembl genes}
#' \description{
#'   This function takes a list of Ensembl gene identifiers which are checked for possible term enrichments.
#'   The fullset is here taken from the Uniprot/Swissprot database
#' }
#' \usage{ enr_enrichment_ensembl(subset, max.genes=5, mapping=NULL, folder=file.path(path.expand("~"),"data")) }
#'
#' \arguments{
#'   \item{subset}{subset of genes, for instance these ones which are upregulated in case of non-human
#'    genes you must submit human Ensembl genes homologous to them}
#'   \item{max.genes}{maximal number of gene identifiers attached into the result table, default: 5}
#'   \item{mapping}{mapping file for any enrichment analysis, must contain a column id and a column gene, if not given GO enrichment is done}
#'   \item{folder}{the data folder where to store intermediate results}
#' }
#' \details{
#'   This function performs an enrichment analysis based on human Ensembl gene identifiers,
#'   if you have non-human genes you should submit the human genes which are homologous to your species genes.
#' }
#' \value{data frame with the following columns 
#' \itemize{
#' \item go_id - the GO or term id 
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
#' genes=read.table(text="ENSG00000102057
#' ENSG00000258659
#' ENSG00000184206
#' ENSG00000184345
#' ENSG00000153684
#' ENSG00000278662
#' ENSG00000183629
#' ENSG00000125520")
#' head(genes[,1])
#' enr$enrichment_ensembl(genes[,1])
#' }
#' 

enr$enrichment_ensembl <- function (subset,max.genes=5,mapping=NULL,folder=file.path(path.expand("~"),"data")) {
   year=as.integer(substr(Sys.Date(),1,4))-1
   dfolder     = folder
   upfile      = uniprot$download("human",folder=dfolder)
   id2gofile   = file.path(dfolder,"humanid2go.tab")
   ensg2gofile = file.path(dfolder,"ensg2go.tab")
   if (is.character(mapping)) {
       mappingfile=mapping
   } else {
       mappingfile = file.path(dfolder,"mapping.tab")
   }
   uniprot$id2go(upfile,id2gofile)   
   uniprot$id2ensg(upfile,ensg2gofile)
   ensg2go=read.table(ensg2gofile,sep="\t")
   colnames(ensg2go)=c("up_id","ensg_id")
   uniprot$mapid2go(id2gofile,ensg2gofile,mappingfile)
   mapping=read.table(mappingfile,sep="\t",quote="",header=TRUE)
   fullset=unique(ensg2go$ensg_id)
   map=mapping[,c('go_id','extern_id')]
   colnames(map)=c("term","gene")
   enrich=enr$enrichment(fullset,subset,map,max.genes=max.genes)
   goutil$new(goutil$download(year,folder=dfolder))
   enrich=cbind(enrich,go_nsp=goutil$getNamespace(enrich$go_id),go_name=goutil$getName(enrich$go_id))
   return(enrich)
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
#' enr$gaf("Apis mellifera", folder=file.path(path.expand("~"),"data"))
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
#'   and converts them into identifiers.
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
#' # convert some terms into unique identifiers
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

#' \name{enr$annotation}
#' \alias{enr_annotation}
#' \alias{enr$annotation}
#' \title{Extract GO annotations for the given genes}
#' \description{
#'   This function takes a vector of gene indentifiers and extracts for them the annotations
#'   from the current GAF file.
#' }
#' \usage{ enr_annotation(x) }
#'
#' \arguments{
#'   \item{x}{vector of gene identifiers}
#' }
#' \details{
#'   This function allows you investigate for a given set of genes what are the available annotations.
#' }
#' \value{data frame with the columns from the annotation file, you can add yourself the GO names as shown in the example}
#' \examples{
#' enr$gaf("Apis mellifera")
#' anno=enr$annotation(c("LOC102653920","LOC725343","DUMMY"))
#' head(anno)
#' fname=goutil$download(2022)
#' 
#' goutil$new(fname)
#' anno=cbind(anno,go_name=goutil$getName(anno$go_id))
#' head(anno)
#' }
#' 

enr$annotation <- function (x) {
    self=enr
    if (is.null(self$.annotation)) {
        stop("Error: You must first initialize your annotations with the enr$gaf method!")
    }
    df=data.frame(gene_id = x) 
    godata=read.table(self$.annotation,comment.char="!",sep="\t",quote="")
    godata=godata[,c(1,2,3,5,6,7,8)]
    colnames(godata)=c("db","id1","id2","go_id","ref","ev_code","nsp")
    if (any(grepl("LOC",x))) {
        godata$id1 = paste("LOC",godata$id1,sep="")
    }
    m1=merge(df,godata,by.x="gene_id",by.y="id1",all.x=TRUE)
    m1=m1[,c("gene_id","db","go_id","ref","ev_code","nsp")]
    m2=merge(df,godata,by.x="gene_id",by.y="id2",all.x=TRUE)
    m2=m2[,c("gene_id","db","go_id","ref","ev_code","nsp")]
    df=unique(rbind(m1,m2))
    return(df)
}

#' \name{enr$lf2dex}
#' \alias{enr_lf2dex}
#' \alias{enr$lf2dex}
#' \title{Create T/F table for a certain log-fold difference}
#' \description{
#'   Create a TRUE/FALSE data frame to indicate if the given two log-fold expressions
#'   against a control for two different conditions are differentially expressed
#' }
#' \usage{ enr_lf2dex(x, threshold.control=1, threshold.groups=0.5) }
#'
#' \arguments{
#'  \item{x}{a two column data frame with log-fold change values for two groups against a control}
#'  \item{threshold.control}{the numerical threshold for the comparison against the control}, 
#'  \item{threshold.groups}{the numerical threshold for the comparison between the groups}against the control, 
#' }
#' \details{
#'   Usually such an analysis can be down using a simple check against a given
#'   log-fold change against a control for instance, if log2 values are compared
#'   data>1, however for the difference set between the two groups, it can happen
#'   that one group has a log-fold change against the control, let's say of 1.3
#'   and the the other group a value of 0.9 against the control. So between each
#'   other they are not differentially expressed. The TRUE/FALSE value in this case
#'   for a  simple check against > 1 would be TRUE (G1) against the control and 
#'   FALSE (G2) against the control.
#' }
#' \value{data frame with two columns with tTRUE or FALSE values for both comparisons}
#' \examples{
#' sdata=data.frame(
#' G1=c(1.2,1.2,1.6,1.2,0.8,-1.6,-1.2,-1.2,-0.1,-0.3,-0.8,0.4,-3.6,-2.6),
#' G2=c(0.9,0.1,1.6,1.8,2.5,-1.4,-0.1,-0.3,-1.2,-1.2,-1.7,1.7,-1.4,-1.4))
#' rownames(sdata)=LETTERS[1:14]
#' sdata
#' ## up-regulated
#' sdata > 1
#' enr$lf2dex(sdata,threshold.control =  1)
#' ## down-regulated
#' sdata < -1
#' enr$lf2dex(sdata,threshold.control = -1)
#' dex=enr$lf2dex(sdata,threshold.control =  1)
#' table(dex[,1],dex[,2])
#' enr$lf2dex(sdata,threshold.control =  1, threshold.groups=0.3)
#' }
#' 


enr$lf2dex <- function (x,threshold.control=1, threshold.groups=0.5) {
    data=x
    # check difference set if it is different between group one and two
    #  as well and not only against the control
    if (threshold.control>0) {
        tab1 = data > threshold.control

        idx1 = tab1[,1] & !(tab1[,2]) & (data[,1]-data[,2] < threshold.groups)
        tab1[idx1,1] = FALSE
        idx2 = tab1[,2] & !(tab1[,1]) & (data[,2]-data[,1] < threshold.groups)
        tab1[idx2,2] = FALSE
    }
    if (threshold.control < 0) {
        print("here")
        tab1 = data < threshold.control
        idx1 = tab1[,1] & !(tab1[,2]) & (data[,1]-data[,2] > threshold.groups)
        tab1[idx1,1] = FALSE
        idx2 = tab1[,2] & !(tab1[,1]) & (data[,2]-data[,1] > threshold.groups)
        tab1[idx2,2] = FALSE
    }
    return(tab1)
}

#' \name{enr$vennplot}
#' \alias{enr_vennplot}
#' \alias{enr$vennplot}
#' \title{Create a two-set venn diagram with the standard color codes}
#' \description{
#'   Use a TRUE/FALSE data frame to create a two set venn diagram plot.
#' }
#' \usage{ enr_vennplot(x, cols=c("#EF536BBB","#61E04FBB","#536BEFBB")) }
#'
#' \arguments{
#'  \item{x}{a two column table with TRUE/FALSE values}
#'  \item{cols}{list of two or three colors with transparency, 
#'              default: c("#EF536BBB","#61E04FBB","#536BEFBB")}
#' }
#' \examples{
#' tab=data.frame(matrix(rnorm(100),ncol=2))
#' colnames(tab)=c('A', 'B')
#' tab=tab>0
#' head(tab)
#' table(tab[,1],tab[,2])
#' enr$vennplot(tab)
#' }
#'

enr$vennplot <- function (x,cols=c("#EF536BBB","#61E04FBB","#536BEFBB")) { 
    venn = function (x,y=NULL,z=NULL,vars=NULL,col=c("#cc888899","#8888cc99","#88cc8899"),cex=1.6,...) {
        circle = function (x,y, radius=1,length=100) {
            theta = seq(0, 2 * pi, length = 100) 
            return(list(x=radius*cos(theta)+x,
                        y=radius*sin(theta)+y))
        }
        venn2D = function (x,col=c("#cc888899","#8888cc99"),cex=1.6,...) {
            if (!is.data.frame(x) & !is.matrix(x)) {
                stop("Error: Not a two column matrix or data frame!")
            }
            if (ncol(x) != 2) {
                stop("Error: Not a two column matrix or data frame!")   
            }
            # reset to useful values, slightly smaller than the defaults:
            # defaults: mai=c(1.02, 0.82, 0.82, 0.42)
            #opar=par(mai=c(1, 0.8, 0.8, 0.4),pty='s')
            # compute circle size
            circ.cex=60*par()$fin[1]/9
            
            plot(c(1,2),c(1,1),xlim=c(0.5,2.5),ylim=c(0.5,2.5),
             pch=19,cex=circ.cex,axes=FALSE,type="n",
             xlab="",ylab="",
             col=col,...)
            polygon(circle(1.2,1.5,radius=0.65),col=col[1],border=col[1])
            polygon(circle(1.8,1.5,radius=0.65),col=col[2],border=col[2])
            text(1.1,2.3,colnames(x)[1],cex=cex)
            text(1.9,2.3,colnames(x)[2],cex=cex)
            # the changes
            if (class(x[,1]) == "logical") {
                is=length(which(x[,1] & x[,2]))
                ls=length(which(x[,1] & !x[,2]))
                rs=length(which(!x[,1] & x[,2]))
                os=length(which(!x[,1] & !x[,2]))
            } else {
                xv=x[,1]
                xv=xv[xv!=""]
                yv=x[,2]
                yv=yv[yv!=""]
                is=length(intersect(xv,yv))
                ls=length(setdiff(xv,yv))
                rs=length(setdiff(yv,xv))
                os=""
            }   
            text(1.5,1.5,is,cex=cex)
            text(0.9,1.5,ls,cex=cex)
            text(2.1,1.5,rs,cex=cex)
            text(1.5,0.7,os,cex=cex)
            #par(opar)
        }                                                   
        if (class(y)[1] != "NULL" & class(z)[1]!="NULL") {
            M=matrix('',ncol=3,nrow=max(c(length(x),length(y),length(z))))
            M[1:length(x),1]=x
            M[1:length(y),2]=y
            M[1:length(z),3]=z        
            colnames(M)=c('x','y','z')
            if (class(vars[1])!="NULL") {
                colnames(M)=vars
            }
            venn(M,col=col,cex=cex,...)
        } else if (class(y)[1] != "NULL") {
            M=matrix('',ncol=2,nrow=max(c(length(x),length(y))))
            M[1:length(x),1]=x
            M[1:length(y),2]=y
            colnames(M)=c('x','y')
            if (class(vars[1])!="NULL") {
                colnames(M)=vars
            }
            venn(M,col=col,cex=cex,...)
        } else if (!is.data.frame(x) & !is.matrix(x)) {
            stop("Error: Not a matrix or data frame!")
        } else if (ncol(x) == 2) {
            venn2D(x,col=col[1:2],cex=cex,...)
        } else if (ncol(x) != 3) {
            stop("Error: Only two or three column matrix or data frame is accepted!")  
        } else if (!class(x[,1]) == "logical") {    
            rnames=unique(c(x[,1],x[,2],x[,3]))
            rnames=rnames[which(rnames!="")]
            M=matrix(FALSE,ncol=3,nrow=length(rnames))
            rownames(M)=rnames
            colnames(M)=colnames(x)
            for (i in 1:3) {
                idx=which(rnames%in%x[,i])
                M[idx,i]=TRUE
            }   
            venn(M,col=col,cex=cex,...)
        } else {  
            #opar=par(mai=c(0.5, 0.4, 0.4, 0.2),pty='s')
            circ.cex=70*par()$fin[1]/9
            plot(c(3.5,5.5,4.5),c(5.5,5.5,3.5),xlim=c(0,9),ylim=c(0,9),
                 pch=19,cex=circ.cex,axes=FALSE,asp=1,
                 xlab="",ylab="", col=col,type="n",...)
            polygon(circle(3.5,5.5,radius=2.3),col=col[1],border=col[1])
            polygon(circle(5.5,5.5,radius=2.3),col=col[2],border=col[2])
            polygon(circle(4.5,3.5,radius=2.3),col=col[3],border=col[3])
            
            text(0.5,7,colnames(x)[1],cex=cex)
            text(8.5,7,colnames(x)[2],cex=cex)
            text(4.5,0.25,colnames(x)[3],cex=cex)
            is=length(which(x[,1] & x[,2] & x[,3]))
            ls=length(which(x[,1] & !x[,2] & !x[,3]))
            rs=length(which(!x[,1] & x[,2] & !x[,3]))
            bs=length(which(!x[,1] & !x[,2] & x[,3]))        
            os=length(which(!x[,1] & !x[,2] & !x[,3]))
            xys=length(which(x[,1] & x[,2] & !x[,3]))
            xzs=length(which(x[,1] & !x[,2] & x[,3]))        
            yzs=length(which(!x[,1] & x[,2] & x[,3]))                
            text(4.5,4.8,is,cex=cex*0.8)
            if (os>0) {
                text(4.5,8.5,os,cex=cex*0.8)
            }
            text(2.2,6,ls,cex=cex*0.8)
            text(6.8,6,rs,cex=cex*0.8)
            text(4.5,2.1,bs,cex=cex*0.8)
            text(4.5,6.6,xys,cex=cex*0.8)    
            text(2.9,4,xzs,cex=cex*0.8)        
            text(6,4,yzs,cex=cex*0.8)            
            #par(opar)
        }   
    }                                                      
    venn(x,col=cols)
}

#' \name{enr$topplot}
#' \alias{enr_topplot}
#' \alias{enr$topplot}
#' \title{Create an color coded plot of expression levels}
#' \description{
#'   Creates a plot were RNA or other expression levels are encoded using more 
#'   or less saturated color schemas.
#' }
#' \usage{ enr_topplot(x, range=c(0,10),legend=TRUE,scale="red",text.start=-0.8) }
#'
#' \arguments{
#'  \item{x}{a two to six column data frame with expression levels, usually on a log2-scale. }
#'  \item{range}{value range used for the scaling of colors}
#'  \item{legend}{Should a color legend be drawn right of the color code rectangles, default: TRUE}
#'  \item{scale}{The color schema used for scaling, possible values are 'red', 'green', 'darkgreen', 'blue' and 'gray' (or 'grey'), default: 'red'} 
#'  \item{text.start}{Where to place te rownames on the x-axis, change this if you have short labels to more positive values, default: -0.8}
#' }
#' \examples{
#' df=data.frame(matrix(runif(104,min=0,max=10),ncol=4))
#' colnames(df)=c("X1","X2", "X3", "X4")
#' rownames(df)=LETTERS[1:26]
#' df=df[order(df[,1]/df[,2]),]
#' par(mfrow=c(1,2),mai=rep(0.1,4))
#' enr$topplot(df,text.start=0.8,scale="darkgreen")
#' enr$topplot(df,text.start=0.8,scale="red")
#' }
#'

enr$topplot <- function (x, range=c(0,10), legend=TRUE, scale="red", text.start=-0.8) {
    mx=ncol(x)
    if (legend) {
        mx=mx+3
    }
    my=nrow(x)
    plot(1,type="n",xlim=c(-1,mx),ylim=c(-2,my+1),xlab="",ylab="",axes=FALSE)
    lastx = 0
    for (xi in 1:ncol(x)) {
        for (yi in 1:nrow(x)) {
            colp=(x[yi,xi]-range[1])/(range[2]-range[1])
            colp=1-colp
            if (scale == "red") {
                col=rgb(1.0,colp,colp)
            } else if (scale == "blue") {
                col=rgb(colp,colp,1.0)
            } else if (scale == "green") {
                col=rgb(colp,1.0,colp)
            } else if (scale == "darkgreen") {
                col=rgb(colp,0.5+colp/2,colp)
            } else if (scale %in% c("gray", "grey")) {
                col=rgb(0.33+colp/1.5,0.33+colp/1.5,0.33+colp/1.5)
            }
            rect(xi*0.5+0.8,my-yi,xi*0.5+1.2,my-yi+1,col=col)
            if (xi == 1) {
                lines(x=c(xi*0.5+0.7,xi*0.5+0.8),y=c(my-yi+0.5,my-yi+0.5))
            }
            lastx = xi*0.5+1.2
        }
    }
    text(text.start,my:1-0.5,rownames(x),pos=4,cex=0.8)
    if (legend) {
        lrange=round(my*0.66)
        lstart=(my-lrange)/2
        #rect(3.45,lstart+1*(lrange/20)-0.05,4.35,lstart+((20+1)*lrange/20)+0.05,border=1,lwd=1)
        for(i in 0:20) {
            colp=1-0.05*i
            if (scale == "red") {
                col=rgb(1.0,colp,colp)
            } else if (scale == "blue") {
                col=rgb(colp,colp,1.0)
            } else if (scale == "green") {
                col=rgb(colp,1.0,colp)
            } else if (scale == "darkgreen") {
                col=rgb(colp,0.5+colp/2,colp)
            } else if (scale == "gray") {
                col=rgb(0.33+colp/1.5,0.33+colp/1.5,0.33+colp/1.5)
            }
            rect(lastx+0.6,lstart+i*(lrange/20),lastx+1.0,lstart+((i+1)*lrange/20),col=col,border=1,lwd=0.2)
            if (i %in% c(0,4,8,12,16,20)) {
                text(lastx+1.1,lstart+(i+0.5)*(lrange/20),sprintf("%0.2f",i*0.05*max(range)),cex=0.8,pos=4)
                lines(x=c(lastx+1.0,lastx+1.1),y=rep(lstart+((i+0.5)*lrange/20),2))

            }
        }
    }
    text(seq(1.5,1.5+(0.5*(ncol(x)-1)),by=0.5),0,paste("|\n",colnames(x),sep=""),pos=1)
}

enr_annotation = enr$annotation
enr_download   = enr$download
enr_enrichment = enr$enrichment
enr_enrichment_ensembl = enr$enrichment_ensembl
enr_gaf        = enr$gaf
enr_name2id    = enr$name2id
enr_new        = enr$new
enr_lf2dex     = enr$lf2dex
enr_symbol2loc = enr$symbol2loc
enr_topplot   = enr$topplot
enr_vennplot   = enr$vennplot

## Private functions


