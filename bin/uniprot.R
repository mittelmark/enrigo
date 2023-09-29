uniprot <- new.env()

uniprot$id2go <- function (infile, outfile) {
    if (grepl("gz$",infile)) {
        fin  = gzfile(infile, "r")
    } else {
        fin  = file(infile, "r")
    }
    fout = file(outfile,'w')
    while(length((line = readLines(fin,n=1)))>0) {
        if (substr(line,1,2)=='ID') {
            id = gsub("^ID\\s+([A-Z0-9_a-z]+).+","\\1",line)
        }  else if (substr(line,1,8) == "DR   GO;") {
            out=gsub("^DR   GO; (GO:[0-9]+); (.):([^;]+); ([A-Z]{2,3}).+","\t\\1\t\\2\t\\3\t\\4",line)
            cat(id,out,"\n",file=fout)            
        }
    }
    close(fin)
    close(fout)
}
uniprot$id2ensg <- function (infile, outfile) {
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
            out=gsub("^DR   Ensembl;.+; (ENSG[0-9]+).+","\t\\1",line)
            cat(id,out,"\n",file=fout)            
        }
    }
    close(fin)
    close(fout)
}
uniprot$mapid2go <- function (upid2gofile,idfile,outfile,evidence="ALL") { 
    id2go=read.table(upid2gofile,sep="\t",quote="");
    colnames(id2go)=c("ID","GO_ID","GO_NAMESPACE","GO_NAME","EVIDENCE")
        
    upid2id=read.table(idfile,sep="\t")
    colnames(upid2id)=c("ID","ID_EXTERN")
    res=merge(id2go,upid2id)
    write.table(res,file=outfile,sep="\t",quote=FALSE)
}
usage  <- function () {
    cat("Usage: uniprot.R --help|--id2ensg|--id2go UNIPROTFILE OUTFILE\n")
    cat("Usage: uniprot.R --merge  UPID2GOFILE ID2EXTERNFILE OUTFILE\n")
}
main <- function (argv) {
    if (length(argv) < 3 | argv[1] %in% c("-h","--help")) {
        usage()
    } else if (argv[1] == "--id2ensg") {
        uniprot$id2ensg(argv[2],argv[3])
    } else if (argv[1] == "--id2go") {
        uniprot$id2go(argv[2],argv[3])
    } else if (argv[1] == "--merge")  {
        if (length(argv) != 4) {
            cat("Error: missing fourth argument (outfile)!")
        } else {
            uniprot$mapid2go(argv[2],argv[3],argv[4])
        }
    }
    
}
if (sys.nframe() == 0L && !interactive()) {
    t1=Sys.time()
    main(commandArgs(trailingOnly=TRUE))
    print(Sys.time()-t1)
}
