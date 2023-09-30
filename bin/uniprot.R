library(enrigo)
usage  <- function () {
    cat("Usage: uniprot.R --help|--i2ensg|--id2go UNIPROTFILE OUTFILE\n")
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
