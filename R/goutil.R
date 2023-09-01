#' \docType{class}  
#' \name{goutil} 
#' \alias{goutil}
#' \alias{goutil-class}
#' \title{ Environment object with functions to work with Gene Ontology terms }
#' \description{
#' This is a collection of useful functions to work with Gene Ontology terms and
#' annotations. Enrichment analysis for gene lists and graphical visualizations are as well provided.#' Below you find the documentation for the following functions 
#' }
#' \section{Methods}{
#' \itemize{
#' \item \code{\link[enrigo:goutil_new]{goutil$new(obofile)}}{initialize the required GO data using the given obo file}
#' \item \code{\link[enrigo:goutil_read.obofile]{goutil$read.obofile(obofile)}}{ reads the given obofile and returns the result}
#' }
#' }
#' \examples{
#' ls(goutil)
#' # get all species
#' fname=goutil$download(2023)
#' print(fname)
#' res=goutil$new(fname)
#' names(res)
#' } 

goutil=new.env()

#' \name{goutil$new}
#' \alias{goutil_new}
#' \alias{goutil$new}
#' \title{ Initializes the goutil object with a new obofile }
#' \description{
#'   This function should be used first to initialize the internal data tables and information.
#' }
#' \usage{ goutil_new(obofile,cache=TRUE) }
#'
#' \arguments{
#'   \item{obofile}{GO obo filename }
#'   \item{cache}{should the results of the parins being cached, default: TRUE }
#' }
#' \details{
#'   This function allows you to intialize your annotation file which contains the mapping between
#'   the GO ids and the gene ids for a given species. For getting a list of available species and for the 
#'   download a internet connection is required.
#' }
#' \value{list with the folloowing components:
#' \itemize{
#'  \item names - with columns id, nsp, and name
#'  \item slims - with columns id, slim
#'  \item obsoletes - with ids which are obsolete
#'  \item tree - with columns child, parent, relation
#' }
#' }
#' \examples{
#' res=goutil$new("2023-01-01-go.obo")
#' sort(names(res))
#' }
#' 

goutil$new <- function (obofile,cache=TRUE) {
    self=goutil
    
    t1=Sys.time()
    
    filename=obofile
    if (!file.exists(filename)) {
        stop("Error: File", filename,"does not exists!")
    }
    self$obofile=obofile
    cache.file=paste(filename,"-data.RDS",sep="")
    if (length(self$godata) == 0 & cache & file.exists(cache.file)) {
        self$gofile=filename
        rds=readRDS(paste(filename,"-data.RDS",sep=""))
        self$godata=rds
        invisible(rds)
    } else if (length(self$godata) != 0 & cache) {
        invisible(self$godata)
    } else {
        res=self$read.obofile(obofile)
        saveRDS(res,file=cache.file)
        self$godata=res
        invisible(res)
    }
}

#' \name{goutil$read.obofile}
#' \alias{goutil_read.obofile}
#' \alias{goutil$read.obofile}
#' \title{ Read the given Gene Ontology Obo file }
#' \description{
#'   This function does the parsing of the obofile, should be used mainly indirectly via `goutil$new`.
#' }
#' \usage{ goutil_read.obofile(obofile) }
#'
#' \arguments{
#'   \item{obofile}{ A GO obo file }
#' }
#' \details{
#'   This function read the actual obofile.
#' }
#' \value{list with the following components:
#' \itemize{
#'  \item names - with columns id, nsp, and name
#'  \item slims - with columns id, slim
#'  \item obsoletes - with ids which are obsolete
#'  \item tree - with columns child, parent, relation
#' }
#' }
#' \examples{
#' res=goutil$read.obofile("2023-01-01-go.obo")
#' sort(names(res))
#' }
#' 

goutil$read.obofile <- function (obofile) {
    if (!file.exists(obofile)) {
        stop("Error: File",obofile,"does not exists!")
    }
    t1=Sys.time()
    tab=read.table(obofile,sep="\t",quote="")
    names=data.frame(id=c(),nsp=c(),name=c())
    tree=data.frame(child=c(),parent=c(),relation=c())
    treel=list()
    sliml=list()
    aidsl=list()
    obs=c()
    idx=grep("^id: GO",tab[,1]);
    for (i in idx) {
        id=tab[i,]
        j=i+1
        while (substr(tab[j,],1,1) != "[") {
            if (grepl("^is_obsolete:",tab[j,1])) {
                obs=c(obs,id)
            } else if (grepl("^is_a:",tab[j,1])) {
                ## rbind is slooooow
                #tree=rbind(tree,data.frame(child=id,parent=tab[j,1],relation="is_a"))
                treel[[j]]=data.frame(child=id,parent=tab[j,1],relation="is_a")
            } else if (grepl("^relationship: part_of",tab[j,1])) {
                treel[[j]]=data.frame(child=id,parent=tab[j,1],relation="part_of")
            } else if (grepl("^subset: goslim",tab[j,1])) {
                sl=gsub("^subset[^_]+_([^ ]+)","\\1",tab[j,1])
                sliml[[j]]=data.frame(id=id,slim=sl)
            } else if (grepl("^alt_id: ",tab[j,1])) { 
                aidsl[[j]]=data.frame(id=id,alt_id=tab[j,1])
            } 
            j=j+1
        }
    }
    ## https://rpubs.com/jimhester/rbind
    tree=do.call(rbind,treel)
    slims=do.call(rbind,sliml)    
    altids=do.call(rbind,aidsl)        
    slims$id=substr(slims$id,5,15)
    tree$child=substr(tree$child,5,15)
    tree$parent=substr(tree$parent,7,16)
    obs=substr(obs,5,15)
    names=data.frame(id=substr(tab[idx,],5,15),nsp=gsub(".+: ","",tab[idx+2,]),name=gsub(".+: ","",tab[idx+1,]))
    print(Sys.time()-t1)
    fin  = file(obofile, "r")
    x=0
    linel=list()
    while(length((line = readLines(fin,n=1)))>0) {
        x=x+1
        if (grepl("^id: GO:",line)) {
            linel[[x]]=list(line,x)
        }
    }
    close(fin)
    lines=do.call(rbind,linel)
    res=list(names=names,obsoletes=obs,
                oboversion=gsub(".+/","",tab[2,1]),
                tree=tree,slims=slims,altids=altids,lines=lines)
    return(res)
}

#' \name{goutil$download}
#' \alias{goutil_download}
#' \alias{goutil$download}
#' \title{ Download for the given year the obofile }
#' \description{
#'   This function is used to download the Obo file for the given year from the gene ontology side.
#' }
#' \usage{ goutil_download(version) }
#'
#' \arguments{
#'   \item{version}{ Either such as 2005 or higher or character string for a version name like "2023-01-01" }
#' }
#' \details{
#'   This function downloads the GO obofile from the Gene Ontolofy website.
#' }
#' \value{filename of the downloaded file}
#' \examples{
#' obofile=goutil$download(2022)
#' }
#' 

goutil$download <- function (version) {
    # versions from 2005-2023 are supperted
    if (is.numeric(version)) {
        if (version == 2022) {
            version=paste(version,"-01-13",sep="")
        } else {
            version=paste(version,"-01-01",sep="")
        }
    }
    if (grepl("^200[0-4]",version)) {
        stop("Only GO obofiles after 2004 can be downloaded!")
    }
    if (!grepl("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$",version)) {
        stop("Error: Invalid version format! Valid versions are 20YY-MM-DD")
    }
    ## http://release.geneontology.org/2005-01-01/ontology/go.obo
    dfile=paste(version,"-go.obo",sep="")
    if (grepl("^2019",version) | grepl("^202",version)) {
        download.file(paste("http://release.geneontology.org/",version,"/ontology/go.obo",sep=""),
                      destfile=dfile)
    } else {
        download.file(paste("http://release.geneontology.org/",version,"/ontology/gene_ontology.obo",sep=""),
                      destfile=dfile)
    }
    return(dfile)
}

goutil_new = goutil$new
goutil_download = goutil$download
goutil_read.obofile = goutil$read.obofile
