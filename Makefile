VERSION := $(shell grep Version: DESCRIPTION | perl -pe 's/.+: //')
PKG     := $(shell basename `pwd`)
build: man/goutil.Rd man/enr.Rd man/uniprot.Rd
	R CMD build .

check: build 
	R CMD check $(PKG)_$(VERSION).tar.gz

man/%.Rd: R/%.R
	Rscript bin/rman.R $<

install: build 
	R CMD INSTALL $(PKG)_$(VERSION).tar.gz
