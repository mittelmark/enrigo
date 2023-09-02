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
#' \item \code{\link[enrigo:goutil_new]{goutil$new}}{initialize the required GO data using the given obo file}
#' \item \code{\link[enrigo:goutil_read.obofile]{goutil$read.obofile}}{ reads the given obofile and returns the result}
#' \item \code{\link[enrigo:goutil_altid2new]{goutil$altid2new}}{convert old alternative GO ids to their new counterpart}
#' \item \code{\link[enrigo:goutil_getChildren]{goutil$getChildren}}{get the child nodes of a given GO id}
#' \item \code{\link[enrigo:goutil_getEntry]{goutil$getEntry}}{get the GO entry in standard text for a given GO id}
#' \item \code{\link[enrigo:goutil_getName]{goutil$getName}}{get the name of (a) given GO ids}
#' \item \code{\link[enrigo:goutil_getNamespace]{goutil$getNamespace}}{get the namespace of (a) given GO ids}
#' \item \code{\link[enrigo:goutil_getOboVersion]{goutil$getOboVersion}}{get the actual version date for the current obo file}
#' \item \code{\link[enrigo:goutil_getParent]{goutil$getParent}}{get the parent node(s) for a given GO id}
#' \item \code{\link[enrigo:goutil_getSlims]{goutil$getSlims}}{get GO-slims, or the GO ids for a given slim}
#' \item \code{\link[enrigo:goutil_getStats]{goutil$getStats}}{statistics for the three main namespaces and the number of obsolete terms}
#' \item \code{\link[enrigo:goutil_getTree]{goutil$getTree}}{reads recursively all parents nodes for a given GO id}
#' \item \code{\link[enrigo:goutil_getTreeMatrix]{goutil$getTreeMatrix}}{adjacency matrix for the given ids in tree}
#' \item \code{\link[enrigo:goutil_isChild]{goutil$isChild}}{check if a given GO id is a direct child of a given parent id}
#' \item \code{\link[enrigo:goutil_isParent]{goutil$isParent}}{check if a given GO id is a direct parent of a given child id}
#' \item \code{\link[enrigo:goutil_kroki]{goutil$kroki}}{create a GO tree graph using the kroki webservie}
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
    lines=as.data.frame(lines)
    colnames(lines)=c('id','line')
    lines$id=substr(lines$id,5,15)
    lines$line=as.numeric(lines$line)
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

#' \name{goutil$altid2new}
#' \alias{goutil_altid2new}
#' \alias{goutil$altid2new}
#' \title{Convert old alternative GO ids to their new counterpart}
#' \description{
#' This function returns the names for a vector of given GO ids. 
#' Before the function can work first a GO obo file must be read in using the function `goutil$new`.
#' }
#' \usage{ goutil_altid2new(goid) }
#'
#' \arguments{
#'   \item{goid}{vector of Gene Ontology ids}
#' }
#' \value{vector of names belonging to the given GO id(s)}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' goutil$altid2new(c("GO:0000001","GO:0061637")) # first should stay unchanged
#' }
#' 

goutil$altid2new <- function (goid) {
    id=goid
    self=goutil
    is_alt = id %in% self$godata$altids$alt_id
    res=c(as.character())
    for (i in 1:length(id)) {
        if (is_alt[i]) {
            res=c(res,as.character(self$godata$altids$id[as.character(self$godata$altids$alt_id)==as.character(id[i])][1]))
        } else {
            res=c(res,as.character(id[i]))
        }
    }
    return(res)
}

#' \name{goutil$getEntry}
#' \alias{goutil_getEntry}
#' \alias{goutil$getEntry}
#' \title{Return the the GO entry as standard text for a given GO id}
#' \description{
#' This function returns the GO entry as a character vector with one line per element.
#' Before the function can work first a GO obo file must be read in using the function `goutil$new`.
#' }
#' \usage{ goutil_getEntry(goid) }
#'
#' \arguments{
#'   \item{goid}{single Gene Ontology ids}
#' }
#' \value{text entry belonging to the given GO id}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' writeLines(unlist(lapply(goutil$getEntry("GO:0000001"),strwrap, 80, exdent=6)))
#' }
#' 

goutil$getEntry = function (goid) {
    self=goutil
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    fin  = file(self$obofile, "r")
    n=self$godata$lines$line[self$godata$lines$id==goid]
    entry=as.character(c())
    readLines(fin,n=n)
    while(length((line = readLines(fin,n=1)))>0) {
        if (grepl('^id:',line)) {
            entry=c(line)
        } else if (grepl('^\\s*$',line)) {
            break
        } else {
            entry=c(entry,line)
        }
        
    }
    close(fin)
    return(entry)
}

#' \name{goutil$getName}
#' \alias{goutil_getName}
#' \alias{goutil$getName}
#' \title{Return the name(s) for given GO id(s)}
#' \description{
#' This function returns a vector of names for all given GO id(s).
#' Before the function can work first a GO obo file must be read in using the function `goutil$new`.
#' }
#' \usage{ goutil_getName(goid) }
#'
#' \arguments{
#'   \item{goid}{vector of Gene Ontology ids}
#' }
#' \details{
#'   This function returns the names for a vector of given GO ids. Before the function
#'   is called you should initialize the goutil environment with the function `goutil$new(obofile)`
#' }
#' \value{vector of names belonging to the given GO id(s)}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' goutil$getName(c("GO:0000001","GO:0000002"))
#' goutil$getName("GO:0003675")
#' }
#' 

goutil$getName <- function (goid) {
    self=goutil
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    if (length(goid)==1) {
        name=self$godata$names$name[self$godata$names$id==goid]
        if (is.null(name)) {
            return(NA)
        } else {
            return(name)
        }
    } else {
        return(unlist(lapply(goid,goutil$getName)))
    }
}

#' \name{goutil$getNamespace}
#' \alias{goutil_getNamespace}
#' \alias{goutil$getNamespace}
#' \title{Return the namespaces(s) for given GO id(s)}
#' \description{
#' This function returns a vector of namespaces for all given GO id(s).
#' Before the function can work first a GO obo file must be read in using the function `goutil$new`.
#' }
#' \usage{ goutil_getNamespace(goid) }
#'
#' \arguments{
#'   \item{goid}{vector of Gene Ontology id(s)}
#' }
#' \details{
#'  This function returns the namespaces with single letter codes for a vector of given GO ids.
#'  Before the function is called you should initialize the goutil environment with the 
#'  function `goutil$new(obofile)`
#' }
#' \value{vector of namespace abbreviations belonging to the given GO id(s)}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' goutil$getNamespace(c("GO:0003675","GO:0008150"))
#' }
#' 

goutil$getNamespace <- function (goid) {
    self=goutil
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    if (length(goid)==1) {
        nsp=self$godata$names$nsp[self$godata$names$id==goid]
        if (is.null(nsp)) {
            return(NA)
        } else {
            return(nsp)
        }
    } else {
        return(unlist(lapply(goid,goutil$getNamespace)))
    }
}
#' \name{goutil$getChildren}
#' \alias{goutil_getChildren}
#' \alias{goutil$getChildren}
#' \title{Return the child nodes of the given GO id}
#' \description{
#'  Return the GO id(s) which have the given go id as a 'is_a' or 'part_of' parent.
#'  Before the function can work first a GO obo file must be read in using the function `goutil$new`.
#' }
#' \usage{ goutil_getChildren(goid) }
#'
#' \arguments{
#'   \item{goid}{single Gene Ontology ids}
#' }
#' \value{vector of child ids having the given GO id as parent}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' goutil$getChildren("GO:0003674")
#' }
#' 

goutil$getChildren <- function (goid) {
    self=goutil
    if (length(self$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    return(self$godata$tree$child[self$godata$tree$parent==goid])
}

#' \name{goutil$getParent}
#' \alias{goutil_getParent}
#' \alias{goutil$getParent}
#' \title{Return the parent GO id(s) for given GO id(s)}
#' \description{
#' This function returns a vector of parent GO id(s) for all given GO id(s).
#' Before the function can work first a GO obo file must be read in using the function `goutil$new`.
#' }
#' \usage{ goutil_getParent(goid,type='all') }
#'
#' \arguments{
#'   \item{goid}{single Gene Ontology ids}
#'   \item{type}{which type of parents should be returned, either 'all', 'is_a' or 'part_of', default: 'all'}
#' }
#' \details{
#'  This function returns the parent GO ids for a vector of given GO ids. This function can be used
#'  for explore the tree of all parent ids for a given GO id.
#'  Before the function is called you should initialize the goutil environment with the 
#'  function `goutil$new(obofile)`
#' }
#' \value{vector of parent ids belonging to the given GO id}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' goutil$getParent("GO:0000001")
#' }
#' 

goutil$getParent <- function (goid,type="all") {
    self=goutil
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    if (type=="all") {
        return(self$godata$tree$parent[self$godata$tree$child == goid])
    } else if (type == "part_of") {
        return(self$godata$tree$parent[self$godata$tree$child == goid & self$godata$tree$relation=="part_of"])
    } else if (type == "is_a") {
        return(self$godata$tree$parent[self$godata$tree$child == goid & self$godata$tree$relation=="is_a"])
    } else {
        stop("Error: Unknown type, must be either 'all', 'part_of' or 'is_a'!")
    }
}

#' \name{goutil$getTree}
#' \alias{goutil_getTree}
#' \alias{goutil$getTree}
#' \title{Return all parent GO id(s) recursively for given GO id(s)}
#' \description{
#' This function returns a vector of all GO ids of parents and their parents
#'   an so on, so recursively.
#' }
#' \usage{ goutil_getTree(goid) }
#'
#' \arguments{
#'   \item{goid}{single Gene Ontology ids}
#' }
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#'   for (id in goutil$getTree("GO:0003676")) { 
#'      cat(id,"\t",goutil$getName(id),"\n") 
#'   }
#' }
#' 

goutil$getTree <- function (goid) {
    self=goutil
    tree=c(goid)
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    for (parent in self$getParent(goid)) {
        tree=c(tree,parent,self$getTree(parent))
    }
    return(unique(tree))
}

#' \name{goutil$getTreeMatrix}
#' \alias{goutil_getTreeMatrix}
#' \alias{goutil$getTreeMatrix}
#' \title{Return an adjacency matrix for the given ids in tree}
#' \description{
#'  This function returns an adjacency matrix for the given list of
#'  GO ids in the form of an adjacency matrix.
#' }
#' \usage{ goutil_getTreeMatrix(tree) }
#'
#' \arguments{
#'   \item{tree}{a list of GO ids, usually create with the goutil$getTee function}
#' }
#' \examples{
#'  goutil$new("2023-01-01-go.obo")
#'  tree=goutil$getTree("GO:0003676")
#'  goutil$getTreeMatrix(tree)
#' }
#'

## Plotting
## https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%7Bids%7D/chart?ids=GO%3A0003676
## > ![See https://www.ebi.ac.uk/QuickGO/term/GO:0003676](img/go0003676-tree.png)


goutil$getTreeMatrix <- function (tree) {
    self=goutil
    M=matrix(0,ncol=length(tree),nrow=length(tree))
    rownames(M)=colnames(M)=tree
    for (i in 1:length(tree)) {
        for (j in 1:length(tree)) {
            if (self$isChild(tree[i],tree[j])) {
                M[tree[j],tree[i]]=1
            }
        }
    }
    return(M)
}

#' \name{goutil$isChild}
#' \alias{goutil_isChild}
#' \alias{goutil$isChild}
#' \title{Checks if the given child is indeed a child node of parent}
#' \description{
#' This function checks for a given child node if it is a 
#' child node of the parent node.
#' }
#' \usage{ goutil_isChild(parent,child) }
#'
#' \arguments{
#'   \item{parent}{tested parent GO id}
#'   \item{child}{tested child GO id}
#' }
#' \value{boolean, TRUE if parent id is a parent node of the child id}
#' \examples{
#'   goutil$new("2023-01-01-go.obo")
#'   goutil$isChild("GO:0003674","GO:0005488") # should be TRUE
#'   goutil$isChild("GO:0003674","GO:0000001") # should be FALSE
#' }
#' 

goutil$isChild <- function (parent,child) {
    return(any(goutil$getChildren(parent) %in% child))
}

#' \name{goutil$isParent}
#' \alias{goutil_isParent}
#' \alias{goutil$isParent}
#' \title{Checks if the given parent node is indeed a parent node of child}
#' \description{
#' This function checks for a given parent node if it is a 
#' parent node of the child node.
#' }
#' \usage{ goutil_isParent(parent,child) }
#'
#' \arguments{
#'   \item{parent}{tested parent GO id}
#'   \item{child}{tested child GO id}
#' }
#' \value{boolean, TRUE if parent id is a parent node of the child id}
#' \examples{
#'   goutil$new("2023-01-01-go.obo")
#'   goutil$isParent("GO:0003674","GO:0005488") # should be TRUE
#'   goutil$isParent("GO:0003674","GO:0000001") # should be FALSE
#' }
#' 

goutil$isParent <- function (parent,child) {
    return(any(goutil$getParent(child) %in% parent))
}


#' \name{goutil$getOboVersion}
#' \alias{goutil_getOboVersion}
#' \alias{goutil$getOboVersion}
#' \title{Returns the actual version date for the current obo file}
#' \description{
#' This function returns a the editing date for the current OBO file in the form YYYY-MM-DD as string.
#' }
#' \usage{ goutil_getOboVersion() }
#'
#' \value{editing date of the Obo file as a string}
#' \examples{
#' obofile=goutil$new("2023-01-01-go.obo")
#' goutil$getOboVersion()
#' }
#' 

goutil$getOboVersion <- function () {
    self=goutil
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    return(self$godata$oboversion)
}

#' \name{goutil$getSlims}
#' \alias{goutil_getSlims}
#' \alias{goutil$getSlims}
#' \title{Return all GO-slims, or the GO ids for a given slim and possibly within the given namespace}
#' \description{
#' This function returns either all existing GO slim names or if
#' if a slim name is given all GO ids, belonging to this slim, if a namespace is given as well the GO ids are checked for the
#' given namespace as well.
#' }
#' \usage{ goutil_getSlims(slim,nsp) }
#'
#' \arguments{
#'   \item{slim}{optional slim name, default: NULL}
#'   \item{nsp}{optional namespace like f (molecular function), c (cellular component), p (biological process), if slim is given and the nsp only GO ids from this namespace will be returned, default: NULL}
#' }
#' \value{vector of slim names or if a slim name is given all GO ids belonging to this slim}
#' \examples{
#' goutil$getSlims()
#' goutil$getSlims(slim="drosophila",nsp="c")
#' for (n in c('c','f','p')) {
#'   print(paste(n,"=",length( goutil$getSlims("drosophila",nsp=n))))
#' }
#' }
#' 
goutil$getSlims <- function (slim=NULL,nsp=NULL) {
    self=goutil
    if (length(self$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    if (is.null(slim)) {
        return(sort(unique(self$godata$slims$slim)))
    } else {
        slims=self$getSlims()
        if (slim %in% slims) {
            ids= self$godata$slims$id[self$godata$slims$slim==slim]
            if (is.null(nsp)) {
                return(ids)
            } else if (nsp %in% c('c','f','p')) {
                set=c()
                tab=self$godata$names[self$godata$names$id %in% ids & self$godata$names$nsp==nsp,]
                return(tab$id)
            } else {
                stop(paste("Error: Unknown namespace",nsp," -  known namespaces are c, f, p!"))
            }
        } else {
            stop(paste("Error: Slim ",slim,"does not exists!"))
        }
    }
}
#' \name{goutil$getStats}
#' \alias{goutil_getStats}
#' \alias{goutil$getStats}
#' \title{Return statistics for the current GO file}
#' \description{
#' This function returns a data frame with a overal summary statistics for the currently loaded GO file.
#' }
#' \usage{ goutil_getStats() }
#'
#' \value{data frame with the columns active and obsolete and the rownames c, f, p for the three namespaces}
#' \examples{
#' goutil$getStats()
#' }
#' 

goutil$getStats <- function () {
    self=goutil
    if (length(goutil$godata)==0) {
        stop("Please read in first GO obo file using goutil$new")
    }
    tab.all=table(goutil$godata$names$nsp)
    tab.obs=table(self$getNamespace(self$godata$obsoletes))
    M=matrix(c(tab.all-tab.obs,tab.obs),ncol=2)
    rownames(M)=names(tab.all)
    colnames(M)=c("valid","obsolete")
    return(M)
}

#' \name{goutil$kroki}
#' \alias{goutil_kroki}
#' \alias{goutil$kroki}
#' \title{Return an URL for a PlantUML graph done using the kroki webservice}
#' \description{
#' This function can be used to embed images for Gene Ontology trees into an HTML page.
#' }
#' \usage{ goutil_kroki(goid) }
#'
#' \arguments{
#'   \item{goid}{single Gene Ontology id}
#' }
#' \value{returns a URL which can be embedded into a Markdown URL for instance}
#' \examples{
#' goutil$kroki("GO:0003676")
#' }
#' 

goutil$kroki <- function (goid) {
    if (!requireNamespace("tcltk")) {
        stop("Funktion goutul$kroki can only be used if package tcltk is installed!")
    }
    tree = goutil$getTree(goid)
    graph2plantuml <- function (g, type="class",letter="C",letter.col="#FFDDEE",dir="forward",labels=c()) {
                                str="
!theme mars        
skinparam MinClassWidth 200
skinparam BackgroundColor #FFFFFF
skinparam defaulttextalignment center
skinparam ClassBackgroundColor #EEFAFF
skinparam ClassHeaderBackgroundColor #EEEEEE

"
        if (length(labels)==0) {
            labels=rep("",ncol(g))
        }
        arr.dir="-->"
        if (dir != "forward") {
            arr.dir = "<--"
        }
        for (i in 1:ncol(g)) {
            lab=paste(strwrap(labels[i],25),collapse="\n")
            str = paste(str,
                        sprintf("class %s << (%s,%s) >> {\n%s\n}\n",colnames(g)[i],letter,letter.col,lab),sep="")
        }
        for (i in 1:(nrow(g)-1)) {
            for (j in (i+1):(nrow(g))) {
                if (g[i,j]==1 & g[j,i] == 1) {
                    str=paste(str,
                              sprintf("\"%s\" -- \"%s\"\n",colnames(g)[i],colnames(g)[j]),sep="")
                } else if (g[i,j]==1 & g[j,i] == 0) {
                    str=paste(str,
                              sprintf("\"%s\" %s \"%s\"\n",colnames(g)[i],arr.dir,colnames(g)[j]),sep="")
                } else if (g[i,j]==0 & g[j,i] == 1) {
                    str=paste(str,
                              sprintf("\"%s\" %s \"%s\"\n",colnames(g)[j],arr.dir,colnames(g)[i]),sep="")
                }                   
            }
        }
        str=paste("@startuml\n",str,"\n@enduml\n")
        return(str)
    }
    plantuml2kroki <- function (text) {
        if (!requireNamespace("tcltk")) {
            stop("Funktion goutil$kroki can only be used if package tcltk is installed!")
        }
        tcltk::.Tcl("
             proc dia2kroki {text} {
                 return [string map {+ - / _ = \"\"}  [binary encode base64 [zlib compress $text]]]
             }
             ")
        url = tcltk::tclvalue(tcltk::tcl("dia2kroki",text))
        url= paste("https://kroki.io/plantuml/png",url,sep="/")
    }

    treem=goutil$getTreeMatrix(tree)
    letter=goutil$getNamespace(tree[1])
    puml=graph2plantuml(treem,letter=toupper(letter),labels=gsub("_"," ",goutil$getName(tree)))
    url=plantuml2kroki(puml)
    return(url)
}

goutil_new = goutil$new
goutil_altid2new = goutil$altid2new
goutil_download = goutil$download
goutil_getChildren = goutil$getChildren
goutil_getEntry = goutil$getEntry
goutil_getName = goutil$getName
goutil_getNamespace = goutil$getNamespace
goutil_getOboVersion = goutil$getOboVersion
goutil_getParent = goutil$getParent
goutil_getSlims = goutil$getSlims
goutil_getStats = goutil$getStats
goutil_getTree = goutil$getTree
goutil_getTreeMatrix = goutil$getTreeMatrix
goutil_isChild = goutil$isChild
goutil_isParent = goutil$isParent
goutil_kroki = goutil$kroki
goutil_read.obofile = goutil$read.obofile
