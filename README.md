# enrigo ![](https://github.com/mittelmark/enrigo/actions/workflows/r.yml/badge.svg)

Enrichment analysis for Gene Ontology terms.

## Functions

Currently the following functions are implemented for the `goutil` object::


- `goutil$new`  - initialize the goutil object
- `goutil$download`  - download a Gene Ontology obo file
- `goutil$read.obofile` - read and parse a Gene Ontology obo file
- `goutil$altid2new` - convert old alternative GO ids to their new counterpart
- `goutil$getChildren` - get the child nodes of a given GO id}
- `goutil$getEntry` - get the GO entry in standard text for a given GO id}
- `goutil$getName` - get the name of (a) given GO ids
- `goutil$getNamespace` - get the namespace of (a) given GO ids
- `goutil$getOboVersion` - get the actual version date for the current obo file
- `goutil$getParent` - get the parent node(s) for a given GO id
- `goutil$getSlims` - get GO-slims, or the GO ids for a given slim
- `goutil$getStats` - statistics for the three main namespaces and the number of obsolete terms}
- `goutil$getTree` - reads recursively all parents nodes for a given GO id
- `goutil$getTreeMatrix` - adjacency matrix for the given ids in tree
- `goutil$isChild` - check if a given GO id is a direct child of a given parent id
- `goutil$isParent` - check if a given GO id is a direct parent of a given child id
- `goutil$kroki` - create a GO tree graph using the kroki webservie
- `goutil$quickGO` - create a GO tree graph using the EBI QuickGO webservie
- `goutil$read.obofile` -  reads the given obofile and returns the result

And  here  the  methods  for the `enr`  object  with  which  you  perform  the
enrichment analysis using GAF files (incomplete):

- `enr$gaf(x)` - initialize a Gene Ontology annotation or query the annotation databases

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

