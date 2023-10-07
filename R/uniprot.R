#' \docType{class}  
#' \name{uniprot} 
#' \alias{uniprot}
#' \alias{uniprot-class}
#' \title{ Environment object with functions to prepare UniProt data for enrichment analysis }
#' \description{
#' The functions within the uniprot environment allow to extract data from UniProt files to do Gene Ontology term enrichment analysis for a given list
#' of genes against a full list of genes from the same organism.
#' }
#' \section{Methods}{
#' \itemize{
#' \item \code{\link[enrigo:uniprot_download]{uniprot$download("taxdivision")}}{download a SwissProt file from the UniProt database from the taxonomic division section}
#' \item \code{\link[enrigo:uniprot_id2go]{uniprot$id2go(infile,outfile)}}{extract mappings between UniProt identifier and GO identifiers}
#' \item \code{\link[enrigo:uniprot_ensg2id]{uniprot$ensg2id2(infile,outfile)}}{extract mappings between UniProt identifiers and Ensembl identifiers}
#' \item \code{\link[enrigo:uniprot_mapid2go]{uniprot$mapid2go(upid2gofile,idfile,outfile,evidence="ALL")}}{create mappings between identifiers and GO identifiers using UniProt identifiers}
#' }
#' }
#' \examples{
#' # download mammmals SwisProt database file
#' file=uniprot$download("human",
#'   folder=file.path(path.expand("~"),"data"))
#' file
#' } 


uniprot = new.env()

#' \name{uniprot$download}
#' \alias{uniprot_download}
#' \alias{uniprot$download}
#' \title{ Download SwissProt UniProt data }
#' \description{
#'   This function downloads SwissProt database files of the UniProt database
#'   to extract from the file GO annotation data.
#' }
#' \usage{ uniprot_download(taxa,timeout=600,folder=NULL) }
#'
#' \arguments{
#'   \item{taxa}{ The taxonomic division, valid values are
#'  "archaea","bacteria","fungi","human","invertebrates","mammals","rodents","vertebrates" or "viruses"}
#'   \item{timeout}{ The allowed timeout in seconds how long the download might need, large downloads might take a long, time, default: 600}
#'   \item{folder}{the data folder where the uniprot data files should be stored, if not given, the current working directory is used, default: NULL}
#' }
#' \details{
#'   This function allows you to prepare annotation file on your own if you do not can ulize the existing
#'   annotation files from the GO consortium. To download the data
#'   a internet connection is required.
#' }
#' \value{local filename, invisible}
#' \examples{
#' # get all species
#' file=uniprot$download("mammals",
#'   folder=file.path(path.expand("~"),"data"))
#' file
#' }
#' 

uniprot$download <- function (taxa,timeout=600,folder=NULL) {
    outfile=sprintf("uniprot_sprot_%s.dat.gz",taxa)
    if (!is.null(folder)) {
        if (!dir.exists(folder)) {
            dir.create(folder,recursive=TRUE)
        }
        outpath=file.path(folder,outfile)
    } else {
        outpath=file.path(getwd(),outfile)
    } 
    if (file.exists(outpath) & substr(file.mtime(outpath),1,7) == substr(Sys.Date(),1,7)) {
        # same month
        invisible(outpath)
    } else {
        to=options()$timeout
        options(timeout=timeout)
        valid_taxa = c("archaea","bacteria","fungi","human","invertebrates","mammals","rodents","vertebrates","viruses")
        if (!(taxa %in% valid_taxa)) {
            stop(paste("Error: Invalid taxa, valid taxa are: '",paste(valid_taxa,collapse="', '"),"'",sep=""))
        }
        download.file(sprintf("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_%s.dat.gz",taxa),outpath)
        options(timeout=to)
        invisible(outpath)
    }
}

#' \name{uniprot$id2go}
#' \alias{uniprot_id2go}
#' \alias{uniprot$id2go}
#' \title{ Extract mappings between uniprot identifier and GO identifiers }
#' \description{
#'   This function creates a mapping file between Uniprot identifiers and GO terms.
#' }
#' \usage{ uniprot_id2go(infile,outfile) }
#'
#' \arguments{
#'   \item{infile}{UniProt file downloaded with `uniprot$download`}
#'   \item{outfile}{file where to store the data}
#' }
#' \details{
#'   This function needs up to several minutes in dependence from the file size.
#' }
#' \value{NULL}
#' \examples{
#' # get all species
#' file=uniprot$download("mammals",
#'    folder=file.path(path.expand("~"),"data"))
#' id2gofile=file.path(path.expand("~"),"data","mammalsid2go.tab")
#' uniprot$id2go(file,id2gofile)
#' head(read.table(id2gofile,sep="\t",quote=""))
#' }
#' 

uniprot$id2go <- function (infile, outfile) {
    if (!(file.exists(outfile) & file.mtime(outfile)>file.mtime(infile))) {
        # only do processing if file does not exists
        # or if infile is newer than outfile
        if (grepl("gz$",infile)) {
            fin  = gzfile(infile, "r")
        } else {
            fin  = file(infile, "r")
        }
        fout = file(outfile,'w')
        cat("uniprot_id\tgo_id\tgo_nsp\tgo_name\tev_code\n",file=fout)
        while(length((line = readLines(fin,n=1)))>0) {
            if (substr(line,1,2)=='ID') {
                id = gsub("^ID\\s+([A-Z0-9_a-z]+).+","\\1",line)
            }  else if (substr(line,1,8) == "DR   GO;") {
                out=gsub("^DR   GO; (GO:[0-9]+); (.):([^;]+); ([A-Z]{2,3}).+","\t\\1\t\\2\t\\3\t\\4",line)
                cat(id,out,"\n",sep="",file=fout)            
            }
        }
        close(fin)
        close(fout)
    }
}

#' \name{uniprot$ensg2id}
#' \alias{uniprot_ensg2id}
#' \alias{uniprot$ensg2id}
#' \title{Extract mappings between uniprot identifier and Ensembl gene indentifier }
#' \description{
#'   This function creates a mapping file between Uniprot identifiers and Ensembl gene identifiers.
#' }
#' \usage{ uniprot_ensg2id(infile,outfile) }
#'
#' \arguments{
#'   \item{infile}{UniProt file downloaded with `uniprot$download`}
#'   \item{outfile}{file where to store the data}
#' }
#' \details{
#'   This function needs up to several minutes in dependence from the file size.
#' }
#' \value{NULL}
#' \examples{
#' # get all species
#' file=uniprot$download("mammals",
#'  folder=file.path(path.expand("~"),"data"))
#' id2ensgfile=file.path(path.expand("~"),"data","id2ensg.tab")
#' uniprot$id2ensg(file,id2ensgfile)
#' head(read.table(id2ensgfile,sep="\t",quote=""))
#' }
#' 

uniprot$id2ensg <- function (infile, outfile) {
    if (!(file.exists(outfile) & file.mtime(outfile)>file.mtime(infile))) {
        # only do processing if file does not exists
        # or if infile is newer than outfile
        
        if (grepl("gz$",infile)) {
            fin  = gzfile(infile, "r")
        } else {
            fin  = file(infile, "r")
        }
        fout = file(outfile,'w')
        while(length((line = readLines(fin,n=1)))>0) {
            if (substr(line,1,2)=='ID') {
                id = gsub("^ID\\s+([A-Z0-9a-z_]+).+","\\1",line)
            }  else if (substr(line,1,13) == "DR   Ensembl;") {
                out=gsub("^DR   Ensembl;.+; (ENS.*G[0-9]+)\\..*","\t\\1",line)
                cat(id,out,"\n",sep="",file=fout)            
            }
        }
        close(fin)
        close(fout)
        tab=read.table(outfile)
        tab=unique(tab)
        colnames(tab)=c("uniprot_id","gene_id")
        write.table(tab,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
    }
}

#' \name{uniprot$mapid2go}
#' \alias{uniprot_mapid2go}
#' \alias{uniprot$mapid2go}
#' \title{Create a mapping between your identifiers and the GO identifiers using identifiers provided by UniProt files}
#' \description{
#'   This function creates a mapping file between your identifiers using
#'   Uniprot identifiers to the GO identifiers provided by the SwissProt database.
#' }
#' \usage{ uniprot_mapid2go(upid2gofile,idfile,outfile,evidence="ALL") }
#'
#' \arguments{
#'   \item{upid2gofile}{file created with `uniprot$id2go`}
#'   \item{idfile}{mapping file for instance created with `uniprot$ensg2go`}
#'   \item{outfile}{outfile name where the mapping is stored, a tabfile}
#'   \item{evidence}{The evidence codes which should be used, currently not used, default: "ALL"}
#' }
#' \value{NULL}
#' \examples{
#' # get all species
#' file=uniprot$download("mammals",
#'   folder=file.path(path.expand("~"),"data"))
#' id2gofile=file.path(path.expand("~"),"data","mammalsid2go.tab")
#' uniprot$id2go(file,id2gofile)
#' id2ensgfile=file.path(path.expand("~"),"data","id2ensg.tab")
#' uniprot$id2ensg(file,id2ensgfile)
#' mappingfile=file.path(path.expand("~"),"data","mapping.tab")
#' uniprot$mapid2go(id2gofile,id2ensgfile,mappingfile)
#' head(read.table(mappingfile,sep="\t",quote=""))
#' }
#' 

uniprot$mapid2go <- function (upid2gofile,idfile,outfile,evidence="ALL") { 
    id2go=read.table(upid2gofile,sep="\t",quote="");
    colnames(id2go)=c("uniprot_id","go_id","go_nsp","go_name","ev_code")
        
    upid2id=read.table(idfile,sep="\t")
    colnames(upid2id)=c("uniprot_id","extern_id")
    upid2id$extern_id=trimws(upid2id$extern_id)
    res=merge(id2go,upid2id)
    write.table(res,file=outfile,sep="\t",quote=FALSE)
}

uniprot_download = uniprot$download
uniprot_id2go    = uniprot$id2go
uniprot_ensg2id  = uniprot$ensg2id
uniprot_mapid2go = uniprot$mapid2go

