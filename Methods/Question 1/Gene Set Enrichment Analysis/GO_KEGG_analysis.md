# GO and KEGG analysis



This file contains quick GO and KEGG analyses using the goana and kegga functions that are included in limma. We will first load the cleaned data and metadata file as well as the "toptables" obtained from the quadratic model in the linear model fitting analysis. The output will not be printed in the md file for better legibility.





## Blood
Perform all analysis for blood first and then repeat it for the other tissue types. First we need to obtain the gene names for each probe. We will do this using the illuminaHumanv4 database in concordance with the platform the samples were run on.

### Get ENTREZ IDs from DE probes in blood

```r
# write code as function to make re-usability easier
add.entrez <- function(TT_Tissue) {
    # obtain probe ids for DE probes in blood:
    probe.ids <- TT_Tissue$probe.id
    # get entrez ID for each probe
    entrez.ids <- unlist(mget(probe.ids, illuminaHumanv4ENTREZID))
    # add to top table
    TT_Tissue$entrez.id <- entrez.ids
    return(TT_Tissue)
}

# apply to blood
TT_Blood <- add.entrez(TT_Blood)
```
We will remove all non-annotated probes:

```r
# number of DE probes in blood:
nrow(TT_Blood)
```

```
## [1] 4103
```

```r
# remove rows that contain non-mapped probes
TT_Blood <- TT_Blood[complete.cases(TT_Blood$entrez.id), ]
# number of probes after removing non-annotated probes
nrow(TT_Blood)
```

```
## [1] 2983
```

```r
length(unique(TT_Blood$entrez.id))
```

```
## [1] 2825
```
After removing non-annotated probes, 2983 of the 4103 probes remain. The 2983 annotated probes correspond to 2825 genes.


```r
# get entrez id for these probes out of the data frame:
blood.DE.entrez <- TT_Blood$entrez.id
```


### GO analysis

```r
blood.go.de <- goana(blood.DE.entrez, trend = TRUE)  #trend=TRUE will adjust for gene abundance
# displaya table with the 20 GO terms with the lowest
# p-value, including only 'biological processes' in the list.
kable(topGO(blood.go.de, ontology = "BP", sort = "DE"), format = "markdown", 
    caption = "Blood GO analysis")
```



|           |Term                                                                |Ont |    N|  DE| P.DE|
|:----------|:-------------------------------------------------------------------|:---|----:|---:|----:|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |   93|  39|    0|
|GO:0045047 |protein targeting to ER                                             |BP  |   99|  39|    0|
|GO:0044033 |multi-organism metabolic process                                    |BP  |  213|  64|    0|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  |  100|  39|    0|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  |  103|  39|    0|
|GO:0019080 |viral gene expression                                               |BP  |  189|  58|    0|
|GO:0019083 |viral transcription                                                 |BP  |  178|  55|    0|
|GO:0044419 |interspecies interaction between organisms                          |BP  |  973| 194|    0|
|GO:0044403 |symbiosis, encompassing mutualism through parasitism                |BP  |  973| 194|    0|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  |  123|  42|    0|
|GO:0016032 |viral process                                                       |BP  |  942| 187|    0|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  |  120|  41|    0|
|GO:0019058 |viral life cycle                                                    |BP  |  444| 103|    0|
|GO:0044764 |multi-organism cellular process                                     |BP  |  947| 187|    0|
|GO:0034655 |nucleobase-containing compound catabolic process                    |BP  |  365|  88|    0|
|GO:0044265 |cellular macromolecule catabolic process                            |BP  |  947| 185|    0|
|GO:0009057 |macromolecule catabolic process                                     |BP  | 1173| 220|    0|
|GO:0019439 |aromatic compound catabolic process                                 |BP  |  406|  94|    0|
|GO:1901361 |organic cyclic compound catabolic process                           |BP  |  429|  97|    0|
|GO:0006401 |RNA catabolic process                                               |BP  |  241|  63|    0|

### KEGG analysis

```r
blood.kegg.de <- kegga(blood.DE.entrez, species = "Hs")
blood.kegg.de <- blood.kegg.de[order(blood.kegg.de$P.DE), ]
blood.kegg.de <- subset(blood.kegg.de, P.DE <= 0.05)
kable(blood.kegg.de, format = "markdown", caption = "Blood KEGG analysis")
```



|              |Pathway            |   N| DE|      P.DE|
|:-------------|:------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome           | 154| 48| 0.0000001|
|path:hsa04742 |Taste transduction |  83| 19| 0.0242738|

Interestingly, the main GO terms associated with the differentially expressed probes in the blood are associated with viral lif cycle (which is unsurprising) but also protein translation/transport and metaoblic processes. No immunology-related GO term is present in the "top20 list". Only two KEGG pathways are significantly associated with the DE probes from blood (p value cutoff = 0.05). 

We will now repeat the same analysis for the other tissues. 

## Colon

```r
# add entrez id tp DE probes in top table:
TT_Colon <- add.entrez(TT_Colon)
# number of DE probes in the colon:
nrow(TT_Colon)
```

```
## [1] 595
```

```r
# remove rows that contain non-mapped probes
TT_Colon <- TT_Colon[complete.cases(TT_Colon$entrez.id), ]
# number of probes after removing non-annotated probes
nrow(TT_Colon)
```

```
## [1] 432
```

```r
length(unique(TT_Colon$entrez.id))
```

```
## [1] 431
```
In the colon, 432 (431 unique genes) of the 595 DE probes are annotated.


```r
# get entrez id for these probes:
colon.DE.entrez <- TT_Colon$entrez.id
# GO analysis
colon.go.de <- goana(colon.DE.entrez, trend = TRUE)
sorted.go.colon <- topGO(colon.go.de, ontology = "BP", sort = "DE")
kable(sorted.go.colon, format = "markdown", caption = "Colon GO analysis")
```



|           |Term                                                                |Ont |   N| DE|      P.DE|
|:----------|:-------------------------------------------------------------------|:---|---:|--:|---------:|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  | 120| 13| 0.0000008|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |  93| 11| 0.0000023|
|GO:0045047 |protein targeting to ER                                             |BP  |  99| 11| 0.0000043|
|GO:0006612 |protein targeting to membrane                                       |BP  | 185| 15| 0.0000046|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  | 100| 11| 0.0000048|
|GO:0019080 |viral gene expression                                               |BP  | 189| 15| 0.0000059|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  | 103| 11| 0.0000064|
|GO:0006413 |translational initiation                                            |BP  | 192| 15| 0.0000072|
|GO:0019083 |viral transcription                                                 |BP  | 178| 14| 0.0000134|
|GO:0006402 |mRNA catabolic process                                              |BP  | 211| 15| 0.0000221|
|GO:0044033 |multi-organism metabolic process                                    |BP  | 213| 15| 0.0000247|
|GO:0006401 |RNA catabolic process                                               |BP  | 241| 16| 0.0000276|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  | 123| 11| 0.0000347|
|GO:0016072 |rRNA metabolic process                                              |BP  | 266| 16| 0.0000897|
|GO:0034655 |nucleobase-containing compound catabolic process                    |BP  | 365| 19| 0.0001386|
|GO:0046700 |heterocycle catabolic process                                       |BP  | 396| 20| 0.0001396|
|GO:0000956 |nuclear-transcribed mRNA catabolic process                          |BP  | 196| 13| 0.0001559|
|GO:0044270 |cellular nitrogen compound catabolic process                        |BP  | 403| 20| 0.0001762|
|GO:0006412 |translation                                                         |BP  | 637| 27| 0.0001925|
|GO:0019439 |aromatic compound catabolic process                                 |BP  | 406| 20| 0.0001943|

```r
# KEGG analysis
colon.kegg.de <- kegga(colon.DE.entrez, species = "Hs")
colon.kegg.de <- colon.kegg.de[order(colon.kegg.de$P.DE), ]
kable(subset(colon.kegg.de, P.DE <= 0.05), format = "markdown", 
    caption = "Colon KEGG analysis")
```



|              |Pathway                                                   |   N| DE|      P.DE|
|:-------------|:---------------------------------------------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome                                                  | 154| 14| 0.0000059|
|path:hsa04961 |Endocrine and other factor-regulated calcium reabsorption |  47|  5| 0.0034024|
|path:hsa04060 |Cytokine-cytokine receptor interaction                    | 270| 12| 0.0147789|
|path:hsa05412 |Arrhythmogenic right ventricular cardiomyopathy (ARVC)    |  72|  5| 0.0201181|
|path:hsa05162 |Measles                                                   | 134|  7| 0.0269759|
|path:hsa05223 |Non-small cell lung cancer                                |  58|  4| 0.0372782|
|path:hsa00790 |Folate biosynthesis                                       |  15|  2| 0.0412551|
|path:hsa05152 |Tuberculosis                                              | 179|  8| 0.0417770|
|path:hsa04721 |Synaptic vesicle cycle                                    |  63|  4| 0.0481936|

## Jejunum

```r
# add entrez id of DE probes in top table:
TT_Jejunum <- add.entrez(TT_Jejunum)
# number of DE probes in the colon:
nrow(TT_Jejunum)
```

```
## [1] 5463
```

```r
# remove rows that contain non-mapped probes
TT_Jejunum <- TT_Jejunum[complete.cases(TT_Jejunum$entrez.id), 
    ]
# number of probes after removing non-annotated probes
nrow(TT_Jejunum)
```

```
## [1] 3977
```

```r
length(unique(TT_Jejunum$entrez.id))
```

```
## [1] 3668
```
After removal of non-annotated probes, 3977 of the 5463 probes remain, corresponding to 3668 genes.


```r
# get entrez id for these probes:
jejunum.DE.entrez <- TT_Jejunum$entrez.id
# GO analysis
jejunum.go.de <- goana(jejunum.DE.entrez, trend = TRUE)
sorted.go.jejunum <- topGO(jejunum.go.de, ontology = "BP", sort = "DE")
kable(sorted.go.jejunum, format = "markdown", caption = "Jejunum GO analysis")
```



|           |Term                                                                |Ont |    N|  DE|     P.DE|
|:----------|:-------------------------------------------------------------------|:---|----:|---:|--------:|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |   93|  38| 0.00e+00|
|GO:0045047 |protein targeting to ER                                             |BP  |   99|  39| 1.00e-07|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  |  103|  40| 1.00e-07|
|GO:0044033 |multi-organism metabolic process                                    |BP  |  213|  66| 4.00e-07|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  |  100|  38| 5.00e-07|
|GO:0019080 |viral gene expression                                               |BP  |  189|  60| 5.00e-07|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  |  120|  43| 6.00e-07|
|GO:0019083 |viral transcription                                                 |BP  |  178|  57| 7.00e-07|
|GO:0016032 |viral process                                                       |BP  |  942| 217| 8.00e-07|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  |  123|  43| 1.20e-06|
|GO:0044764 |multi-organism cellular process                                     |BP  |  947| 217| 1.20e-06|
|GO:0044419 |interspecies interaction between organisms                          |BP  |  973| 222| 1.30e-06|
|GO:0044403 |symbiosis, encompassing mutualism through parasitism                |BP  |  973| 222| 1.30e-06|
|GO:0033036 |macromolecule localization                                          |BP  | 2749| 553| 3.00e-06|
|GO:0034655 |nucleobase-containing compound catabolic process                    |BP  |  365|  95| 8.50e-06|
|GO:0072657 |protein localization to membrane                                    |BP  |  478| 117| 1.86e-05|
|GO:0006401 |RNA catabolic process                                               |BP  |  241|  67| 1.90e-05|
|GO:0008104 |protein localization                                                |BP  | 2379| 477| 2.38e-05|
|GO:0006612 |protein targeting to membrane                                       |BP  |  185|  54| 2.77e-05|
|GO:0006605 |protein targeting                                                   |BP  |  677| 156| 2.88e-05|

```r
# KEGG analysis
jejunum.kegg.de <- kegga(jejunum.DE.entrez, species = "Hs")
jejunum.kegg.de <- jejunum.kegg.de[order(jejunum.kegg.de$P.DE), 
    ]
kable(subset(jejunum.kegg.de, P.DE <= 0.05), format = "markdown", 
    caption = "Jejunum KEGG analysis")
```



|              |Pathway                                                |   N| DE|      P.DE|
|:-------------|:------------------------------------------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome                                               | 154| 47| 0.0000677|
|path:hsa03450 |Non-homologous end-joining                             |  13|  8| 0.0005382|
|path:hsa04120 |Ubiquitin mediated proteolysis                         | 137| 39| 0.0012481|
|path:hsa05200 |Pathways in cancer                                     | 395| 92| 0.0025672|
|path:hsa00010 |Glycolysis / Gluconeogenesis                           |  67| 21| 0.0047717|
|path:hsa05414 |Dilated cardiomyopathy                                 |  90| 26| 0.0061732|
|path:hsa05412 |Arrhythmogenic right ventricular cardiomyopathy (ARVC) |  72| 21| 0.0116608|
|path:hsa05410 |Hypertrophic cardiomyopathy (HCM)                      |  83| 23| 0.0161393|
|path:hsa05145 |Toxoplasmosis                                          | 113| 29| 0.0219331|
|path:hsa04010 |MAPK signaling pathway                                 | 255| 58| 0.0236230|
|path:hsa00561 |Glycerolipid metabolism                                |  59| 17| 0.0246485|
|path:hsa04914 |Progesterone-mediated oocyte maturation                |  96| 25| 0.0268190|
|path:hsa00604 |Glycosphingolipid biosynthesis - ganglio series        |  15|  6| 0.0364536|
|path:hsa04261 |Adrenergic signaling in cardiomyocytes                 | 144| 34| 0.0445028|
|path:hsa04721 |Synaptic vesicle cycle                                 |  63| 17| 0.0448203|
|path:hsa04975 |Fat digestion and absorption                           |  41| 12| 0.0483716|

## Tonsil

```r
# add entrez id of DE probes in top table:
TT_Tonsil <- add.entrez(TT_Tonsil)
# number of DE probes in the colon:
nrow(TT_Tonsil)
```

```
## [1] 316
```

```r
# remove rows that contain non-mapped probes
TT_Tonsil <- TT_Tonsil[complete.cases(TT_Tonsil$entrez.id), ]
# number of probes after removing non-annotated probes
nrow(TT_Tonsil)
```

```
## [1] 233
```

```r
length(unique(TT_Tonsil$entrez.id))
```

```
## [1] 232
```
After removal of non-annotated probes, 233 (232 unique genes) of the 316 probes remain.


```r
# get entrez id for these probes:
tonsil.DE.entrez <- TT_Tonsil$entrez.id
# GO analysis
tonsil.go.de <- goana(tonsil.DE.entrez, trend = TRUE)
sorted.go.tonsil <- topGO(tonsil.go.de, ontology = "BP", sort = "DE")
kable(sorted.go.tonsil, format = "markdown", caption = "Tonsil GO analysis")
```



|           |Term                                                                |Ont |    N| DE|      P.DE|
|:----------|:-------------------------------------------------------------------|:---|----:|--:|---------:|
|GO:0006612 |protein targeting to membrane                                       |BP  |  185| 11| 0.0000053|
|GO:0006605 |protein targeting                                                   |BP  |  677| 21| 0.0000134|
|GO:0006886 |intracellular protein transport                                     |BP  | 1011| 24| 0.0002254|
|GO:0006413 |translational initiation                                            |BP  |  192|  9| 0.0002367|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  |  120|  7| 0.0003172|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |   93|  6| 0.0004998|
|GO:0045047 |protein targeting to ER                                             |BP  |   99|  6| 0.0006967|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  |  100|  6| 0.0007347|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  |  103|  6| 0.0008582|
|GO:0033036 |macromolecule localization                                          |BP  | 2749| 46| 0.0012189|
|GO:1902582 |single-organism intracellular transport                             |BP  |  699| 17| 0.0014778|
|GO:0044765 |single-organism transport                                           |BP  | 2700| 45| 0.0015144|
|GO:0090150 |establishment of protein localization to membrane                   |BP  |  354| 11| 0.0016000|
|GO:0072594 |establishment of protein localization to organelle                  |BP  |  643| 16| 0.0016108|
|GO:0071702 |organic substance transport                                         |BP  | 2584| 43| 0.0020329|
|GO:0072661 |protein targeting to plasma membrane                                |BP  |   24|  3| 0.0020915|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  |  123|  6| 0.0021381|
|GO:0034638 |phosphatidylcholine catabolic process                               |BP  |    7|  2| 0.0023262|
|GO:0034613 |cellular protein localization                                       |BP  | 1559| 29| 0.0025459|
|GO:0070727 |cellular macromolecule localization                                 |BP  | 1570| 29| 0.0028201|

```r
# KEGG analysis
tonsil.kegg.de <- kegga(tonsil.DE.entrez, species = "Hs")
tonsil.kegg.de <- tonsil.kegg.de[order(tonsil.kegg.de$P.DE), 
    ]
kable(subset(tonsil.kegg.de, P.DE <= 0.05), format = "markdown", 
    caption = "Tonsil KEGG analysis")
```



|              |Pathway              |   N| DE|      P.DE|
|:-------------|:--------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome             | 154|  7| 0.0009376|
|path:hsa00340 |Histidine metabolism |  24|  2| 0.0246020|

## Axillary Lymph Node

```r
# add entrez id of DE probes in top table:
TT_ALN <- add.entrez(TT_ALN)
# number of DE probes in the colon:
nrow(TT_ALN)
```

```
## [1] 2561
```

```r
# remove rows that contain non-mapped probes
TT_ALN <- TT_ALN[complete.cases(TT_ALN$entrez.id), ]
# number of probes after removing non-annotated probes
nrow(TT_ALN)
```

```
## [1] 1831
```

```r
length(unique(TT_ALN$entrez.id))
```

```
## [1] 1777
```
After removal of non-annotated probes, 1831 of the 2561 probes remain, resulting in 1777 uniqe entrez IDs.


```r
# get entrez id for these probes:
aln.DE.entrez <- TT_ALN$entrez.id
# GO analysis
aln.go.de <- goana(aln.DE.entrez, trend = TRUE)
sorted.go.aln <- topGO(aln.go.de, ontology = "BP", sort = "DE")
kable(sorted.go.aln, format = "markdown", caption = "Axillary Lymph Node GO analysis")
```



|           |Term                                                                |Ont |     N|   DE|      P.DE|
|:----------|:-------------------------------------------------------------------|:---|-----:|----:|---------:|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |    93|   23| 0.0000011|
|GO:0045047 |protein targeting to ER                                             |BP  |    99|   23| 0.0000036|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  |   100|   23| 0.0000043|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  |   103|   23| 0.0000073|
|GO:0072657 |protein localization to membrane                                    |BP  |   478|   65| 0.0000318|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  |   120|   24| 0.0000331|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  |   123|   24| 0.0000504|
|GO:0006612 |protein targeting to membrane                                       |BP  |   185|   31| 0.0000960|
|GO:0090150 |establishment of protein localization to membrane                   |BP  |   354|   50| 0.0000979|
|GO:0033365 |protein localization to organelle                                   |BP  |   878|  103| 0.0001092|
|GO:0019083 |viral transcription                                                 |BP  |   178|   30| 0.0001110|
|GO:0044033 |multi-organism metabolic process                                    |BP  |   213|   34| 0.0001216|
|GO:0019080 |viral gene expression                                               |BP  |   189|   31| 0.0001438|
|GO:0034397 |telomere localization                                               |BP  |    12|    6| 0.0001765|
|GO:0006413 |translational initiation                                            |BP  |   192|   31| 0.0001926|
|GO:0006402 |mRNA catabolic process                                              |BP  |   211|   33| 0.0002247|
|GO:0044699 |single-organism process                                             |BP  | 12914| 1115| 0.0002299|
|GO:0072594 |establishment of protein localization to organelle                  |BP  |   643|   78| 0.0002650|
|GO:0009057 |macromolecule catabolic process                                     |BP  |  1173|  129| 0.0002691|
|GO:0006401 |RNA catabolic process                                               |BP  |   241|   36| 0.0003008|

```r
# KEGG analysis
aln.kegg.de <- kegga(aln.DE.entrez, species = "Hs")
aln.kegg.de <- aln.kegg.de[order(aln.kegg.de$P.DE), ]
kable(subset(aln.kegg.de, P.DE <= 0.05), format = "markdown", 
    caption = "Axillary Lymph Node KEGG analysis")
```



|              |Pathway                                                    |   N| DE|      P.DE|
|:-------------|:----------------------------------------------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome                                                   | 154| 31| 0.0000026|
|path:hsa00601 |Glycosphingolipid biosynthesis - lacto and neolacto series |  27|  7| 0.0054005|
|path:hsa00780 |Biotin metabolism                                          |   3|  2| 0.0195653|
|path:hsa00591 |Linoleic acid metabolism                                   |  29|  6| 0.0293765|
|path:hsa00603 |Glycosphingolipid biosynthesis - globo and isoglobo series |  15|  4| 0.0308210|
|path:hsa05412 |Arrhythmogenic right ventricular cardiomyopathy (ARVC)     |  72| 11| 0.0343962|
|path:hsa04722 |Neurotrophin signaling pathway                             | 119| 16| 0.0368322|
|path:hsa03060 |Protein export                                             |  23|  5| 0.0374414|

## Mesenteric Lymph Node

```r
# add entrez id of DE probes in top table:
TT_MLN <- add.entrez(TT_MLN)
# number of DE probes in the colon:
nrow(TT_MLN)
```

```
## [1] 2063
```

```r
# remove rows that contain non-mapped probes
TT_MLN <- TT_MLN[complete.cases(TT_MLN$entrez.id), ]
# number of probes after removing non-annotated probes
nrow(TT_MLN)
```

```
## [1] 1489
```

```r
length(unique(TT_MLN$entrez.id))
```

```
## [1] 1450
```
After removal of non-annotated probes, 1489 (1450 unique genes) of the 2063 probes remain.


```r
# get entrez id for these probes:
mln.DE.entrez <- TT_MLN$entrez.id
# GO analysis
mln.go.de <- goana(mln.DE.entrez, trend = TRUE)
sorted.go.mln <- topGO(mln.go.de, ontology = "BP", sort = "DE")
kable(sorted.go.mln, format = "markdown", caption = "Mesenteric Lymph Node GO analysis")
```



|           |Term                                                                |Ont |    N|  DE|     P.DE|
|:----------|:-------------------------------------------------------------------|:---|----:|---:|--------:|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  |  100|  25| 0.00e+00|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |   93|  24| 0.00e+00|
|GO:0045047 |protein targeting to ER                                             |BP  |   99|  24| 0.00e+00|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  |  103|  24| 1.00e-07|
|GO:0019058 |viral life cycle                                                    |BP  |  444|  58| 6.00e-07|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  |  120|  24| 1.00e-06|
|GO:0009057 |macromolecule catabolic process                                     |BP  | 1173| 120| 1.30e-06|
|GO:0006413 |translational initiation                                            |BP  |  192|  32| 1.30e-06|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  |  123|  24| 1.60e-06|
|GO:0006612 |protein targeting to membrane                                       |BP  |  185|  31| 1.70e-06|
|GO:1901566 |organonitrogen compound biosynthetic process                        |BP  | 1348| 133| 2.30e-06|
|GO:0006518 |peptide metabolic process                                           |BP  |  804|  86| 8.10e-06|
|GO:0043603 |cellular amide metabolic process                                    |BP  |  940|  97| 9.70e-06|
|GO:1901575 |organic substance catabolic process                                 |BP  | 1838| 168| 1.02e-05|
|GO:0090150 |establishment of protein localization to membrane                   |BP  |  354|  46| 1.03e-05|
|GO:0034655 |nucleobase-containing compound catabolic process                    |BP  |  365|  47| 1.06e-05|
|GO:0072657 |protein localization to membrane                                    |BP  |  478|  57| 1.34e-05|
|GO:0019083 |viral transcription                                                 |BP  |  178|  28| 1.85e-05|
|GO:0019439 |aromatic compound catabolic process                                 |BP  |  406|  50| 1.91e-05|
|GO:0046700 |heterocycle catabolic process                                       |BP  |  396|  49| 2.04e-05|

```r
# KEGG analysis
mln.kegg.de <- kegga(mln.DE.entrez, species = "Hs")
mln.kegg.de <- mln.kegg.de[order(mln.kegg.de$P.DE), ]
kable(subset(mln.kegg.de, P.DE <= 0.05), format = "markdown", 
    caption = "Mesenteric Lymph Node KEGG analysis")
```



|              |Pathway                                                |   N| DE|      P.DE|
|:-------------|:------------------------------------------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome                                               | 154| 28| 0.0000033|
|path:hsa05412 |Arrhythmogenic right ventricular cardiomyopathy (ARVC) |  72| 11| 0.0124911|
|path:hsa04975 |Fat digestion and absorption                           |  41|  7| 0.0244558|
|path:hsa04725 |Cholinergic synapse                                    | 112| 14| 0.0277756|
|path:hsa04120 |Ubiquitin mediated proteolysis                         | 137| 16| 0.0343578|
|path:hsa04340 |Hedgehog signaling pathway                             |  47|  7| 0.0474815|

## Genital Pelvic Lymph Node

```r
# add entrez id of DE probes in top table:
TT_GLN <- add.entrez(TT_GLN)
# number of DE probes in the colon:
nrow(TT_GLN)
```

```
## [1] 1886
```

```r
# remove rows that contain non-mapped probes
TT_GLN <- TT_GLN[complete.cases(TT_GLN$entrez.id), ]
# number of probes after removing non-annotated probes
nrow(TT_GLN)
```

```
## [1] 1366
```

```r
length(unique(TT_GLN$entrez.id))
```

```
## [1] 1338
```
After removal of non-annotated probes, 1366 of the 1886 probes remain, corresponding to 1338 unique genes.


```r
# get entrez id for these probes:
gln.DE.entrez <- TT_GLN$entrez.id
# GO analysis
gln.go.de <- goana(gln.DE.entrez, trend = TRUE)
sorted.go.gln <- topGO(gln.go.de, ontology = "BP", sort = "DE")
kable(sorted.go.gln, format = "markdown", caption = "Genital-Pelvic Lymph Node GO analysis")
```



|           |Term                                                                |Ont |    N|  DE|     P.DE|
|:----------|:-------------------------------------------------------------------|:---|----:|---:|--------:|
|GO:0006413 |translational initiation                                            |BP  |  192|  34| 0.00e+00|
|GO:0006614 |SRP-dependent cotranslational protein targeting to membrane         |BP  |   93|  20| 8.00e-07|
|GO:0006401 |RNA catabolic process                                               |BP  |  241|  35| 2.10e-06|
|GO:0045047 |protein targeting to ER                                             |BP  |   99|  20| 2.20e-06|
|GO:0006613 |cotranslational protein targeting to membrane                       |BP  |  100|  20| 2.60e-06|
|GO:0009057 |macromolecule catabolic process                                     |BP  | 1173| 111| 3.80e-06|
|GO:0000184 |nuclear-transcribed mRNA catabolic process, nonsense-mediated decay |BP  |  120|  22| 3.90e-06|
|GO:0072599 |establishment of protein localization to endoplasmic reticulum      |BP  |  103|  20| 4.20e-06|
|GO:0070972 |protein localization to endoplasmic reticulum                       |BP  |  123|  22| 5.90e-06|
|GO:0043043 |peptide biosynthetic process                                        |BP  |  663|  70| 7.80e-06|
|GO:0006518 |peptide metabolic process                                           |BP  |  804|  81| 8.90e-06|
|GO:0043604 |amide biosynthetic process                                          |BP  |  732|  75| 1.07e-05|
|GO:0006412 |translation                                                         |BP  |  637|  67| 1.37e-05|
|GO:0006402 |mRNA catabolic process                                              |BP  |  211|  30| 1.69e-05|
|GO:0022613 |ribonucleoprotein complex biogenesis                                |BP  |  462|  52| 2.09e-05|
|GO:0043603 |cellular amide metabolic process                                    |BP  |  940|  90| 2.10e-05|
|GO:0042981 |regulation of apoptotic process                                     |BP  | 1358| 121| 2.39e-05|
|GO:0000956 |nuclear-transcribed mRNA catabolic process                          |BP  |  196|  28| 2.93e-05|
|GO:0043067 |regulation of programmed cell death                                 |BP  | 1370| 121| 3.54e-05|
|GO:0034655 |nucleobase-containing compound catabolic process                    |BP  |  365|  43| 3.66e-05|

```r
# KEGG analysis
gln.kegg.de <- kegga(gln.DE.entrez, species = "Hs")
gln.kegg.de <- gln.kegg.de[order(gln.kegg.de$P.DE), ]
kable(subset(gln.kegg.de, P.DE <= 0.05), format = "markdown", 
    caption = "Genital-Pelvic Lymph Node KEGG analysis")
```



|              |Pathway                                                |   N| DE|      P.DE|
|:-------------|:------------------------------------------------------|---:|--:|---------:|
|path:hsa03010 |Ribosome                                               | 154| 26| 0.0000068|
|path:hsa05412 |Arrhythmogenic right ventricular cardiomyopathy (ARVC) |  72| 11| 0.0069103|
|path:hsa05410 |Hypertrophic cardiomyopathy (HCM)                      |  83| 11| 0.0193213|
|path:hsa04919 |Thyroid hormone signaling pathway                      | 116| 14| 0.0194257|
|path:hsa04722 |Neurotrophin signaling pathway                         | 119| 14| 0.0237351|
|path:hsa04120 |Ubiquitin mediated proteolysis                         | 137| 15| 0.0348825|
|path:hsa00591 |Linoleic acid metabolism                               |  29|  5| 0.0385087|

## Conclusion
In conclusion, the main GO biological processes associated with the DE genes after SIV infection are related to protein translation and metabolic processes as well as viral life cycle. Viral replication requires translation using the host's ribosome, increasing the energy demand within infected host cells. The results suggests rapid replication of the virus in host cells throughout the entire body. 

Why the DE probes from several tissues are significantly enriched for Arrhythmogenic right ventricular cardiomyopathy is not obvious and would require more analysis of the genes that drive the enrichment. 

Interestingly, the only tissue where an immune-related KEGG pathway is significantly enriched is the colon ("Cytokine-cytokine receptor interaction" pathway). In general, there are significantly less immune cells present in the colon compared to all other tissues analyzed. The enrichment for cytokine-cytokine receptor interaction in the colon might indicate infiltration of immune cells (in response to infection of colon cells). This would have to be confirmed on a functional level. Alternatively, colon cells might upregulate cytokine receptors in response to circulating cytokines.

It is striking, that the KEGG pathway "ribosome" appears as significantly enriched (with very low p values) across all tissue types. This has potential implications for the significance of the microarray results. The transcriptomic results suggest major changes in protein translation upon SIV infection. If true, the translation of different mRNA transcripts might be affected to a different extent (e.g. increased translation of mRNAs with certain 3'UTR structures but not others). In how far the transcriptomic changes translate into difference on the protein level is unclear in every study but will be complicated by changes in the ribosomal machinery. Therefore, functional studies to strengthen the findings from this microarray analysis are esstential.
