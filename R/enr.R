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

enr_gaf = enr$gaf
enr_new = enr$new
