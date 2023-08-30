VERSION := $(shell grep Version: DESCRIPTION | perl -pe 's/.+: //')
PKG     := $(shell basename `pwd`)
build:
	R CMD build .

check: build man/add.Rd
	echo $(VERSION)
	R CMD check $(PKG)_$(VERSION).tar.gz

man/%.Rd: R/%.R
	Rscript bin/rman.R $<
