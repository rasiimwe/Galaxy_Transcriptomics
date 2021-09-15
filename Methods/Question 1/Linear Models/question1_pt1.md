# Question 1: How does gene expression change in various tissues after SIV infection?

## Fitting Linear Models and Identifying DEG


The tissues selected for analysis were the jejunum, blood, tonsil, axillary lymph node, mesenteric lymph node, genital pelvic lymph node, and colon because each had at least 3 biological replicates for all of the time points. From model validataion analysis of our exploratiry fitting we find that models including a quadratic term are a better fit to the data. Therefore we need to fit tissue speciifc linear models that include the quadratic of time. 




```r
## Load cleaned data and metadata
DATA <- read.table("../team_SIV-in-Rhesus-Monkeys/Data/Processed Data/DATA.txt", 
    header = TRUE, row.names = 1)
colnames(DATA) <- gsub("X", "", colnames(DATA))
MetaData <- read.table("../team_SIV-in-Rhesus-Monkeys/Data/Processed Data/MetaData_cleaned.txt", 
    header = TRUE, row.names = 1)
## Subsetting the metadata file for each tissue of interest
## and sorting each subset metadata by day.
MetaData_Jejunum <- subset(MetaData, MetaData$tissue == "Jejunum")
MetaData_Jejunum <- MetaData_Jejunum[with(MetaData_Jejunum, order(days)), 
    ]
```
Now that we have our subset data frame for the Jejunum we can fit our linear model.

```r
## Design a matrix considering the variable time and include a
## quadratic of time to fit a non-linear model:
Design_Jejunum <- model.matrix(~days + I(days^2), data = MetaData_Jejunum)
## Subset the total data with all of the Jejunum samples to
## match the matrix:
DATA_Jejunum <- DATA[, c(row.names(MetaData_Jejunum))]
## Fit the data to the desgn matrix:
Model_Jejunum <- eBayes(lmFit(DATA_Jejunum, Design_Jejunum))
## Get the toptable of the model:
TT_Jejunum <- topTable(Model_Jejunum, coef = 2:3, n = Inf)
## Because the filtering function gets rid of the row names we
## need to set the rownames as a seperate column to not loose
## the gene.id
TT_Jejunum$probe.id <- row.names(TT_Jejunum)
## Subset the toptable by FDR and p_value thresholds
TT_Jejunum <- (TT_Jejunum %>% filter(adj.P.Val <= 0.05))
```
How many probes are found to be differentially expressed?

```r
nrow(TT_Jejunum)
```

```
## [1] 5463
```

## The same methodology is used for the other tissues: 
### Blood

```r
MetaData_Blood <- subset(MetaData, MetaData$tissue == "Blood")
MetaData_Blood <- MetaData_Blood[with(MetaData_Blood, order(days)), 
    ]
DATA_Blood <- DATA[, c(row.names(MetaData_Blood))]
Design_Blood <- model.matrix(~days + I(days^2), data = MetaData_Blood)
Model_Blood <- eBayes(lmFit(DATA_Blood, Design_Blood))
TT_Blood <- topTable(Model_Blood, coef = 2:3, n = Inf)
TT_Blood$probe.id <- row.names(TT_Blood)
TT_Blood <- (TT_Blood %>% filter(adj.P.Val <= 0.05))
nrow(TT_Blood)
```

```
## [1] 4103
```

### Tonsil 

```r
MetaData_Tonsil <- subset(MetaData, MetaData$tissue == "Tonsil")
MetaData_Tonsil <- MetaData_Tonsil[with(MetaData_Tonsil, order(days)), 
    ]
DATA_Tonsil <- DATA[, c(row.names(MetaData_Tonsil))]
Design_Tonsil <- model.matrix(~days + I(days^2), data = MetaData_Tonsil)
Model_Tonsil <- eBayes(lmFit(DATA_Tonsil, Design_Tonsil))
TT_Tonsil <- topTable(Model_Tonsil, coef = 2:3, n = Inf)
TT_Tonsil$probe.id <- row.names(TT_Tonsil)
TT_Tonsil <- (TT_Tonsil %>% filter(adj.P.Val <= 0.05))
nrow(TT_Tonsil)
```

```
## [1] 316
```

### Axillary Lymph Node 

```r
MetaData_ALN <- subset(MetaData, MetaData$tissue == "axillary_LN")
MetaData_ALN <- MetaData_ALN[with(MetaData_ALN, order(days)), 
    ]
DATA_ALN <- DATA[, c(row.names(MetaData_ALN))]
Design_ALN <- model.matrix(~days + I(days^2), data = MetaData_ALN)
Model_ALN <- eBayes(lmFit(DATA_ALN, Design_ALN))
TT_ALN <- topTable(Model_ALN, coef = 2:3, n = Inf)
TT_ALN$probe.id <- row.names(TT_ALN)
TT_ALN <- (TT_ALN %>% filter(adj.P.Val <= 0.05))
nrow(TT_ALN)
```

```
## [1] 2561
```

### Mesenteric Lymph Node 

```r
MetaData_MLN <- subset(MetaData, MetaData$tissue == "mesenteric_LN")
MetaData_MLN <- MetaData_MLN[with(MetaData_MLN, order(days)), 
    ]
DATA_MLN <- DATA[, c(row.names(MetaData_MLN))]
Design_MLN <- model.matrix(~days + I(days^2), data = MetaData_MLN)
Model_MLN <- eBayes(lmFit(DATA_MLN, Design_MLN))
TT_MLN <- topTable(Model_MLN, coef = 2:3, n = Inf)
TT_MLN$probe.id <- row.names(TT_MLN)
TT_MLN <- (TT_MLN %>% filter(adj.P.Val <= 0.05))
nrow(TT_MLN)
```

```
## [1] 2063
```

### Genital Pelvic Lymph Node 

```r
MetaData_GLN <- subset(MetaData, MetaData$tissue == "genital_pelvic_LN")
MetaData_GLN <- MetaData_GLN[with(MetaData_GLN, order(days)), 
    ]
DATA_GLN <- DATA[, c(row.names(MetaData_GLN))]
Design_GLN <- model.matrix(~days + I(days^2), data = MetaData_GLN)
Model_GLN <- eBayes(lmFit(DATA_GLN, Design_GLN))
TT_GLN <- topTable(Model_GLN, coef = 2:3, n = Inf)
TT_GLN$probe.id <- row.names(TT_GLN)
TT_GLN <- (TT_GLN %>% filter(adj.P.Val <= 0.05))
nrow(TT_GLN)
```

```
## [1] 1886
```

### Colon 

```r
MetaData_Colon <- subset(MetaData, MetaData$tissue == "Colon")
MetaData_Colon <- MetaData_Colon[with(MetaData_Colon, order(days)), 
    ]
DATA_Colon <- DATA[, c(row.names(MetaData_Colon))]
Design_Colon <- model.matrix(~days + I(days^2), data = MetaData_Colon)
Model_Colon <- eBayes(lmFit(DATA_Colon, Design_Colon))
TT_Colon <- topTable(Model_Colon, coef = 2:3, n = Inf)
TT_Colon$probe.id <- row.names(TT_Colon)
TT_Colon <- (TT_Colon %>% filter(adj.P.Val <= 0.05))
nrow(TT_Colon)
```

```
## [1] 595
```

# Conclusion: 

There are differences in the number of differentially expressed genes between tissues through time. The table below summarizes the number of hits in each tissue of interest: 

| Model | Number of hits at FDR <= 0.05 |
| ------------- | ------------- |
| Jejunum | 5463 | 
| Blood | 4103 | 
| Tonsil | 316 | 
| Auxillary Lymph Node | 2561 | 
| Mesenteric Lymph Node | 2063 | 
| Genital Pelvic Lymph Node | 1886 | 
| Colon | 595 | 
