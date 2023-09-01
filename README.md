# enrigo ![](https://github.com/mittelmark/enrigo/actions/workflows/r.yml/badge.svg)

Enrichment analysis for Gene Ontology terms.

## Functions

Currently the following functions are implemented:


- `goutil$new(obofile)`  - initialize the goutil object
- `goutil$download(version)`  - download a Gene Ontology obo file
- `goutil$read.obofile(obofile)` - read and parse a Gene Ontology obo file
- `gaf$new(x)` - initialize a Gene Ontology annotation

## Installation 

Currently the package is a WIP, it is not yet advisable to install it.

You could install it with the library remotes like this:

```r
> install.libarary(remotes)
> library(remotes)
> remotes::install_github("https://github.com/mittelmark/enrigo")
```

## Author and Copyright

Authors: 

- Detlef Groth, University of Potsdam, Germany
- Nahal Zaman, University of Potsdam, Germany
- Merhshad Shobeiri, University of Potsdam, Germany

License: MIT License see the file [LICENSE](LICENSE) for details.

## Bug reporting

In case of bugs and suggestions, use the [issues](https://github.com/mittelmark/enrigo/issues) link on top.

[issues](../../issues)
