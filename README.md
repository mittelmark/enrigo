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

And  here  the  methods  for the `enr`  object  with  which  you  perform  the
enrichment analysis using GAF files (incomplete):

- `enr$gaf` - initialize a Gene Ontology annotation or query the annotation
  databases
- `enr$new` - initialize the `enr` object with the download annotation file
- `enr$enrichment` - perform the enrichment analysis
- `enr$symbol2loc` - convert gene symbols into NCBI gene identifiers

## Installation 

You could install it with the library remotes like this:

```r
> install.libarary(remotes)
> library(remotes)
> remotes::install_github("https://github.com/mittelmark/enrigo")
```

## Enrichment analysis

After you installed the package you might have a look at the package  vignette
like this:

```r
vignette("tutorial",package="enrigo")
```

Here in short the steps to perform an enrichment analysis:

* load the enrigo package and other packages you might need
* create a list of all genes for which you have data (fullset)
* create a subset of this list of genes which are you genes of  interest,  for
  instance up- or down-regulated genes (subset)
* initialize the `goutil` object with a specific Gene Ontology Obo file version
* initialize the `enr` object with you species name
* call the `enr$enrichment` function with the fullset and the subset vectors

Here is some R code assuming you have gene expression data for `Apis melifera`
and your gene identifiers are in the variables `fullset` and `subset`

```r
library(enrigo)
fname=goutil$download(2022)
goutil$new(fname)
enr$gaf('Apis mellifera')
enr$new()
enr$enrichment(fullset,subset)
```

The  `enr$enrichment`  function will create a table with columns  representing
the GO  terms,  the  p-values  for over or under  representation,  statistical
measures  like effect sizes,  adjusted  p-values etc. For more details look at
the package vignette.

## Supported Species

The following  species are annotated with gene ontology  annotation files (you
can obtain the list by using the `enr$gaf` function): 


| Species                      | Annotation File                                                                                                                            |                                                     
|:-----------------------------|:----------------------------------------------------------------------------------------------------------------------------|                                                     
|Acromyrmex echinatior         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Acromyrmex_echinatior_HGD_go_annotation.gaf.gz        |                                                     
|Apis cerana                   |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Apis_cerana_HGD_go_annotation.gaf.gz                  |                                                     
|Apis dorsata                  |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Apis_dorsata_HGD_go_annotation.gaf.gz                 |                                                     
|Apis florea                   |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Apis_florea_HGD_go_annotation.gaf.gz                  |                                                     
|Apis mellifera                |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Apis_mellifera_HGD_go_annotation.gaf.gz               |                                                     
|Arabidopsis thaliana          |http://current.geneontology.org/annotations/tair.gaf.gz                                                                      |                                                     
|Athalia rosae                 |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Athalia_rosae_HGD_go_annotation.gaf.gz                |                                                     
|Atta cephalotes               |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Atta_cephalotes_HGD_go_annotation.gaf.gz              |
|Atta colombica                |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Atta_colombica_HGD_go_annotation.gaf.gz               |
|Belonocnema treatae           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Belonocnema_treatae_HGD_go_annotation.gaf.gz          |
|Bombus bifarius               |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Bombus_bifarius_HGD_go_annotation.gaf.gz              |
|Bombus impatiens              |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Bombus_impatiens_HGD_go_annotation.gaf.gz             |
|Bombus terrestris             |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Bombus_terrestris_HGD_go_annotation.gaf.gz            |
|Bombus vancouverensis         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Bombus_vancouverensis_HGD_go_annotation.gaf.gz        |
|Bombus vosnesenskii           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Bombus_vosnesenskii_HGD_go_annotation.gaf.gz          |
|Caenorhabditis elegans        |http://current.geneontology.org/annotations/wb.gaf.gz                                                                        |                                                     
|Camponotus floridanus         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Camponotus_floridanus_HGD_go_annotation.gaf.gz        |
|Cardiocondyla obscurior       |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Cardiocondyla_obscurior_HGD_go_annotation.gaf.gz      |
|Cephus cinctus                |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Cephus_cinctus_HGD_go_annotation.gaf.gz               |
|Ceratina calcarata            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Ceratina_calcarata_HGD_go_annotation.gaf.gz           |
|Ceratosolen solmsi            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Ceratosolen_solmsi_HGD_go_annotation.gaf.gz           |
|Chelonus insularis            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Chelonus_insularis_HGD_go_annotation.gaf.gz           |
|Copidosoma floridanum         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Copidosoma_floridanum_HGD_go_annotation.gaf.gz        |
|Cyphomyrmex costatus          |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Cyphomyrmex_costatus_HGD_go_annotation.gaf.gz         |
|Danio rerio                   |http://current.geneontology.org/annotations/zfin.gaf.gz                                                                      |                                                     
|Diachasma alloeum             |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Diachasma_alloeum_HGD_go_annotation.gaf.gz            |
|Dictyostelium discoideum      |http://current.geneontology.org/annotations/dictybase.gaf.gz                                                                 |                                                     
|Dinoponera quadriceps         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Dinoponera_quadriceps_HGD_go_annotation.gaf.gz        |
|Drosophila melanogaster       |http://current.geneontology.org/annotations/fb.gaf.gz                                                                        |                                                     
|Dufourea novaeangliae         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Dufourea_novaeangliae_HGD_go_annotation.gaf.gz        |
|Escherichia coli              |http://current.geneontology.org/annotations/ecocyc.gaf.gz                                                                    |                                                     
|Eufriesea mexicana            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Eufriesea_mexicana_HGD_go_annotation.gaf.gz           |
|Fopius arisanus               |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Fopius_arisanus_HGD_go_annotation.gaf.gz              |
|Formica exsecta               |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Formica_exsecta_HGD_go_annotation.gaf.gz              |
|Habropoda laboriosa           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Habropoda_laboriosa_HGD_go_annotation.gaf.gz          |
|Harpegnathos saltator         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Harpegnathos_saltator_HGD_go_annotation.gaf.gz        |
|Homo sapiens                  |http://current.geneontology.org/annotations/cgd.gaf.gz                                                                       |                                                     
|Lasioglossum albipes          |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Lasioglossum_albipes_HGD_go_annotation.gaf.gz         |
|Leishmania major              |http://current.geneontology.org/annotations/xenbase.gaf.gz                                                                   |                                                     
|Linepithema humile            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Linepithema_humile_HGD_go_annotation.gaf.gz           |
|Megachile rotundata           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Megachile_rotundata_HGD_go_annotation.gaf.gz          |
|Megalopta genalis             |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Megalopta_genalis_HGD_go_annotation.gaf.gz            |
|Melipona quadrifasciata       |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Melipona_quadrifasciata_HGD_go_annotation.gaf.gz      |
|Microplitis demolitor         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Microplitis_demolitor_HGD_go_annotation.gaf.gz        |
|Monomorium pharaonis          |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Monomorium_pharaonis_HGD_go_annotation.gaf.gz         |
|Mus musculus                  |http://current.geneontology.org/annotations/sgn.gaf.gz                                                                       |                                                     
|Nasonia vitripennis           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Nasonia_vitripennis_HGD_go_annotation.gaf.gz          |
|Neodiprion lecontei           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Neodiprion_lecontei_HGD_go_annotation.gaf.gz          |
|Nomia melanderi               |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Nomia_melanderi_HGD_go_annotation.gaf.gz              |
|Nylanderia fulva              |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Nylanderia_fulva_HGD_go_annotation.gaf.gz             |
|Odontomachus brunneus         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Odontomachus_brunneus_HGD_go_annotation.gaf.gz        |
|Ooceraea biroi                |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Ooceraea_biroi_HGD_go_annotation.gaf.gz               |
|Orussus abietinus             |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Orussus_abietinus_HGD_go_annotation.gaf.gz            |
|Osmia bicornis                |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Osmia_bicornis_HGD_go_annotation.gaf.gz               |
|Osmia lignaria                |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Osmia_lignaria_HGD_go_annotation.gaf.gz               |
|Plasmodium falciparum         |http://current.geneontology.org/annotations/genedb_pfalciparum.gaf.gz                                                        |                                                     
|Pogonomyrmex barbatus         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Pogonomyrmex_barbatus_HGD_go_annotation.gaf.gz        |
|Polistes canadensis           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Polistes_canadensis_HGD_go_annotation.gaf.gz          |
|Polistes dominula             |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Polistes_dominula_HGD_go_annotation.gaf.gz            |
|Pseudomonas aeruginosa        |http://current.geneontology.org/annotations/pseudocap.gaf.gz                                                                 |                                                     
|Pseudomyrmex gracilis         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Pseudomyrmex_gracilis_HGD_go_annotation.gaf.gz        |
|Rattus norvegicus             |http://current.geneontology.org/annotations/rgd.gaf.gz                                                                       |                                                     
|Saccharomyces cerevisiae      |http://current.geneontology.org/annotations/sgd.gaf.gz                                                                       |                                                     
|Schizosaccharomyces japonicus |http://current.geneontology.org/annotations/reactome.gaf.gz                                                                  |                                                     
|Schizosaccharomyces pombe     |http://current.geneontology.org/annotations/pombase.gaf.gz                                                                   |                                                     
|Solenopsis invicta            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Solenopsis_invicta_HGD_go_annotation.gaf.gz           |
|Temnothorax curvispinosus     |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Temnothorax_curvispinosus_HGD_go_annotation.gaf.gz    |
|Trachymyrmex cornetzi         |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Trachymyrmex_cornetzi_HGD_go_annotation.gaf.gz        |
|Trachymyrmex septentrionalis  |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Trachymyrmex_septentrionalis_HGD_go_annotation.gaf.gz |
|Trachymyrmex zeteki           |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Trachymyrmex_zeteki_HGD_go_annotation.gaf.gz          |
|Trichogramma pretiosum        |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Trichogramma_pretiosum_HGD_go_annotation.gaf.gz       |
|Trypanosoma brucei            |http://current.geneontology.org/annotations/genedb_tbrucei.gaf.gz                                                            |                                                     
|Vollenhovia emeryi            |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Vollenhovia_emeryi_HGD_go_annotation.gaf.gz           |
|Wasmannia auropunctata        |https://data.elsiklab.missouri.edu/downloads/hgd/HGD-GO-Annotation/gaf/Wasmannia_auropunctata_HGD_go_annotation.gaf.gz       |

If you know other repositories of gene ontology  annotation files you can tell
me, please use the Bug reporting link at the bottom.


## Author and Copyright

Authors: 

- Detlef Groth, University of Potsdam, Germany
- Nahal Zaman, University of Potsdam, Germany
- Merhshad Shobeiri, University of Potsdam, Germany

License: MIT License see the file [LICENSE](LICENSE) for details.

## Bug reporting

In      case      of      bugs       and       suggestions,       use      the
[issues](https://github.com/mittelmark/enrigo/issues) link on top.



