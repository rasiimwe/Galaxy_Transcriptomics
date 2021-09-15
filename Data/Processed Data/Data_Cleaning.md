# Data Cleaning




## Data Import:

We are working with a publicly available data set. After download from the GEO website (GSE80013), we will now inspect and clean the data.

Download supplemental files including raw data:

```r
filePaths = getGEOSuppFiles("GSE80013")
```

```
## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80013/suppl/
```

```
## OK
```

```r
# load and inspect the non-normalized matrix:
rawDATA <- read.table("GSE80013/GSE80013_non-normalized_matrix.txt.gz", 
    header = TRUE, row.names = 1)
kable(rawDATA[1:4, 1:5], format = "markdown", digits = 2)
```



|             | X8292044049_A| X8292044049_B| X8292044049_C| X8292044049_D| X8292044049_E|
|:------------|-------------:|-------------:|-------------:|-------------:|-------------:|
|ILMN_1802380 |           312|           458|           124|           477|           147|
|ILMN_1893287 |            88|            92|            85|            86|            81|
|ILMN_3238331 |            87|            80|            80|            90|            88|
|ILMN_1736104 |            85|            90|            86|            97|            79|

```r
# get dimensions of the data frame:
dim(rawDATA)
```

```
## [1] 47231   697
```

The data frame contains 697 columns (samples) and 47231 rows (probes) as expected. How many NA values are in the dataframe?

```r
sum(is.na(rawDATA))
```

```
## [1] 1034
```

This number is relatively low for such a big data frame. We will remove all rows from the dataset that contain missing values.

```r
rawDATA <- na.omit(rawDATA)
```

Read the metadata file:

```r
MetaData <- read.table("../team_SIV-in-Rhesus-Monkeys/Data/Raw Data/Metadata.txt", 
    header = TRUE)
# sanity check: display excerpt of the meta data:
kable(head(MetaData), format = "markdown")
```



|             | subject_ID|tissue      | days|
|:------------|----------:|:-----------|----:|
|8292044049_A | 8292044049|Inguinal_LN |    7|
|8292044049_B | 8292044049|Liver       |    7|
|8292044049_C | 8292044049|Jejunum     |    7|
|8292044049_D | 8292044049|Blood       |    7|
|8292044049_E | 8292044049|Thymus      |    7|
|8292044049_F | 8292044049|Tonsil      |    7|

```r
# change the column names of the data file to match the
# metadata
colnames(rawDATA) <- gsub("X", "", colnames(rawDATA))
```

## Selection of relevant tissue samples:
We are only interested in the seven tissues for which there is at least three replicates at each time point (Jejunum, Blood, Tonsil, Colon as well as the axiallary, mesenteric and genital-pelvic Lymph Node). We will select these samples and drop all other samples from the MetaData and the rawDATA files.

```r
# make character string containing the tissues of interest
tissues_of_interest <- c("Jejunum", "Blood", "Tonsil", "axillary_LN", 
    "mesenteric_LN", "genital_pelvic_LN", "Colon")
# subset the MetaData keeping only the samples from tissues
# of interest
MetaData <- MetaData[MetaData$tissue %in% tissues_of_interest, 
    ]
MetaData <- droplevels(MetaData)
# extract the samples names of the selected samples
samples_of_interest <- row.names(MetaData)
# subset the data keeping only the samples of interest
rawDATA <- rawDATA[, samples_of_interest]
# get dimensions of the data frame
dim(rawDATA)
```

```
## [1] 47028   231
```

After selection of the tissues of interest, 231 samples remain in out data set.

Sanity check: check for log2 intensity distribution between the samples:

```r
lograwDATA <- log2(rawDATA)
boxplot(lograwDATA, range = 0, ylab = "log2 intensity", xaxt = "n", 
    main = "sample log2 intensity distribution")
```

<img src ="https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/unnamed-chunk-5-1.png">

## Sample to sample correlation:

```r
# correlate samples and store in variable c:
c <- cor(lograwDATA)
```

We will look at the distribution of inter-sample variability:

```r
hist(c, sub = paste("Mean=", format(mean(c[upper.tri(c)]), digits = 3)), 
    main = "inter-sample variability")
```
<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/inter-sample%20variability-1.png">

```r
min(c)
```

```
## [1] 0.5689715
```

The tail to the left suggests that there might be some outliers present within the data. The lowest correlation value is 0.57. 

We will visualize the inter-sample correlation in a heatmap, ordering the samples by tissue and day post infection:

```r
# arrange the dataset by tissue and within each tissue by
# days post infection
tissue_day <- MetaData[order(MetaData$tissue, MetaData$days), 
    ]
tissue_day_names <- rownames(tissue_day)
# create a heatmap displaying the correlation between the
# samples:
aheatmap(c[tissue_day_names, tissue_day_names], Rowv = NA, Colv = NA, 
    cellwidth = 1.2, cellheight = 1.2, annRow = list(Tissue = tissue_day$tissue, 
        Day = as.factor(tissue_day$days)))
```

<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/correlation%20heatmap-1.png">

It looks like there is high correlation within each tissue independent of the time point. Also, there is high correlation between the three lymph node tissues and blood, which is expected as theses tissues types are physiologically very similar. Based on the heatmap, some outliers seem to be present. To assess whether this is due to biological variability (tissue differences), we will create a heatmap for tonsil, jejunum and colon each. We will create a heatmap containing axillary, mesenteric and genital-pelvic lymphnode and blood. We will group these tissue samples in one heatmap as there seems to be little variation between them which is physiologically reasonable. The color scale will be the same for all heatmaps to facilitate comparison.

### Colon

```r
colon <- MetaData[MetaData$tissue == "Colon", ]  #subset meta data by tissue of interest
colon <- colon[order(colon$days), ]  #oder by days
colon_names <- rownames(colon)  #retrieve sample names
aheatmap(c[colon_names, colon_names], Rowv = NA, Colv = NA, cellwidth = 5.5, 
    cellheight = 5.5, annRow = list(Day = as.factor(colon$days)), 
    breaks = 0.8)  #make heatmap of the subsetted correlation matrix
```

<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/colon%20correlation-1.png">

Sample 8909358261_A seems to be an outlier. 

```r
# get summary statistic for every sample
kable(t(apply(c[colon_names, colon_names], 2, summary)), format = "markdown")
```



|             |   Min.| 1st Qu.| Median|   Mean| 3rd Qu.| Max.|
|:------------|------:|-------:|------:|------:|-------:|----:|
|9982865061_K | 0.7540|  0.9134| 0.9316| 0.9201|  0.9421|    1|
|9982865068_F | 0.7448|  0.9097| 0.9286| 0.9180|  0.9386|    1|
|9982865068_H | 0.7587|  0.9195| 0.9397| 0.9272|  0.9494|    1|
|9372962047_B | 0.7456|  0.9510| 0.9631| 0.9537|  0.9728|    1|
|9377358071_H | 0.7853|  0.9146| 0.9219| 0.9170|  0.9369|    1|
|9377361031_C | 0.7408|  0.9391| 0.9486| 0.9411|  0.9570|    1|
|9377361031_D | 0.7550|  0.9460| 0.9647| 0.9537|  0.9734|    1|
|8909358010_L | 0.7495|  0.9333| 0.9623| 0.9485|  0.9728|    1|
|9287078038_I | 0.7479|  0.9422| 0.9649| 0.9516|  0.9745|    1|
|9287078041_C | 0.7543|  0.9475| 0.9675| 0.9541|  0.9761|    1|
|9287078041_H | 0.7020|  0.9279| 0.9442| 0.9340|  0.9560|    1|
|9287078048_K | 0.7476|  0.9488| 0.9605| 0.9482|  0.9682|    1|
|9377358032_K | 0.7283|  0.9496| 0.9579| 0.9496|  0.9669|    1|
|9377358053_E | 0.6856|  0.9193| 0.9363| 0.9255|  0.9550|    1|
|9377358053_L | 0.6602|  0.8789| 0.8985| 0.8961|  0.9272|    1|
|9377361036_L | 0.6634|  0.8821| 0.9004| 0.8998|  0.9349|    1|
|8381688093_G | 0.7039|  0.9331| 0.9463| 0.9358|  0.9552|    1|
|8909358001_I | 0.7547|  0.9339| 0.9535| 0.9431|  0.9704|    1|
|8909358010_D | 0.7653|  0.9331| 0.9557| 0.9436|  0.9673|    1|
|8909358015_H | 0.7254|  0.9300| 0.9557| 0.9443|  0.9665|    1|
|8909358026_L | 0.7332|  0.9387| 0.9607| 0.9466|  0.9671|    1|
|8909358125_F | 0.7616|  0.9408| 0.9614| 0.9503|  0.9753|    1|
|8909358275_A | 0.7568|  0.9430| 0.9702| 0.9545|  0.9775|    1|
|8981245020_J | 0.7278|  0.9433| 0.9524| 0.9446|  0.9591|    1|
|8981245021_D | 0.7442|  0.9363| 0.9549| 0.9414|  0.9579|    1|
|8981245042_G | 0.6656|  0.8990| 0.9150| 0.9102|  0.9433|    1|
|9021212001_E | 0.7365|  0.9301| 0.9432| 0.9344|  0.9521|    1|
|9031313015_E | 0.7489|  0.9427| 0.9612| 0.9491|  0.9732|    1|
|9031313032_G | 0.7486|  0.9485| 0.9620| 0.9500|  0.9697|    1|
|8909358015_L | 0.7220|  0.9511| 0.9596| 0.9491|  0.9678|    1|
|8909358261_A | 0.6602|  0.7278| 0.7456| 0.7427|  0.7547|    1|
|8909358294_C | 0.7356|  0.9397| 0.9618| 0.9479|  0.9705|    1|
|9021212001_L | 0.7561|  0.9421| 0.9522| 0.9404|  0.9628|    1|

This supports that sample 8909358261_A is an outlier. It has a mean correlation of 0.74 with the other colon samples. All other colon samples have a mean correlation of 0.89 or higher. We therefore remove this sample from our dataset.

Remove sample 8909358261_A from data and metadata files:

```r
rawDATA <- rawDATA[, -as.numeric(colnames(rawDATA) == "8909358261_A")]
MetaData <- MetaData[-as.numeric(rownames(MetaData) == "8909358261_A"), 
    ]
```

### Tonsil

```r
tonsil <- MetaData[MetaData$tissue == "Tonsil", ]
tonsil <- tonsil[order(tonsil$days), ]
tonsil_names <- rownames(tonsil)
aheatmap(c[tonsil_names, tonsil_names], Rowv = NA, Colv = NA, 
    cellwidth = 9, cellheight = 9, annRow = list(Day = as.factor(tonsil$days)), 
    breaks = 0.8)
```

<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/tonsil%20correlation-1.png">

Sample 9982865061_I seems to be an outlier and sample 9377358065_K might be an outlier as well.


```r
# get summary statistic for every sample
kable(t(apply(c[tonsil_names, tonsil_names], 2, summary)), format = "markdown")
```



|             |   Min.| 1st Qu.| Median|   Mean| 3rd Qu.| Max.|
|:------------|------:|-------:|------:|------:|-------:|----:|
|9980097039_A | 0.8140|  0.9300| 0.9426| 0.9365|  0.9507|    1|
|9980097039_E | 0.7995|  0.9218| 0.9455| 0.9362|  0.9549|    1|
|9982865061_I | 0.6839|  0.7737| 0.8028| 0.8067|  0.8373|    1|
|9982865068_I | 0.7893|  0.8758| 0.9045| 0.8978|  0.9188|    1|
|9372962066_B | 0.8432|  0.9390| 0.9525| 0.9462|  0.9630|    1|
|9372962066_H | 0.7219|  0.9318| 0.9524| 0.9367|  0.9610|    1|
|9377358063_A | 0.8400|  0.9347| 0.9595| 0.9476|  0.9639|    1|
|9377361008_D | 0.8316|  0.9515| 0.9621| 0.9536|  0.9706|    1|
|8909358005_G | 0.7811|  0.9488| 0.9646| 0.9524|  0.9737|    1|
|9287078026_C | 0.7844|  0.9484| 0.9643| 0.9509|  0.9709|    1|
|9287078026_J | 0.7516|  0.9426| 0.9520| 0.9417|  0.9630|    1|
|9287078036_G | 0.7197|  0.9278| 0.9546| 0.9370|  0.9630|    1|
|9287078048_D | 0.7910|  0.9467| 0.9661| 0.9521|  0.9717|    1|
|9372535050_D | 0.8139|  0.9561| 0.9624| 0.9551|  0.9746|    1|
|9377358019_G | 0.8030|  0.9501| 0.9657| 0.9537|  0.9726|    1|
|9377358065_K | 0.6839|  0.8281| 0.8405| 0.8370|  0.8486|    1|
|8292044049_F | 0.7175|  0.9294| 0.9499| 0.9339|  0.9598|    1|
|8909358001_B | 0.8025|  0.9544| 0.9658| 0.9551|  0.9727|    1|
|8909358006_G | 0.7507|  0.9463| 0.9611| 0.9479|  0.9709|    1|
|8909358014_H | 0.8144|  0.9430| 0.9603| 0.9496|  0.9649|    1|
|8909358017_H | 0.7828|  0.9415| 0.9575| 0.9460|  0.9672|    1|
|8909358017_I | 0.8015|  0.9556| 0.9684| 0.9560|  0.9728|    1|
|8981245016_K | 0.6853|  0.9125| 0.9378| 0.9207|  0.9453|    1|
|8981245022_I | 0.7998|  0.9543| 0.9659| 0.9548|  0.9731|    1|
|8981245037_G | 0.7974|  0.9549| 0.9650| 0.9545|  0.9712|    1|
|9031313012_K | 0.8086|  0.9529| 0.9592| 0.9520|  0.9709|    1|
|9031313013_I | 0.8306|  0.9511| 0.9624| 0.9534|  0.9676|    1|
|9031313023_E | 0.7854|  0.9148| 0.9359| 0.9299|  0.9512|    1|
|8909358015_J | 0.8465|  0.9502| 0.9636| 0.9529|  0.9698|    1|
|8909358294_D | 0.8303|  0.9241| 0.9382| 0.9326|  0.9473|    1|
|9031313012_J | 0.7474|  0.9351| 0.9498| 0.9363|  0.9579|    1|
|9031313013_G | 0.8156|  0.9398| 0.9563| 0.9458|  0.9617|    1|

Both have a mean correlation with the other tonsil samples that is much lower than the mean values for the other tonsil samples, so we will remove these outliers from the dataset.

```r
rawDATA <- rawDATA[, -as.numeric(colnames(rawDATA) == "9982865061_I")]
MetaData <- MetaData[-as.numeric(rownames(MetaData) == "9982865061_I"), 
    ]

rawDATA <- rawDATA[, -as.numeric(colnames(rawDATA) == "9377358065_K")]
MetaData <- MetaData[-as.numeric(rownames(MetaData) == "9377358065_K"), 
    ]
```

### Jejunum

```r
jejunum <- MetaData[MetaData$tissue == "Jejunum", ]
jejunum <- jejunum[order(jejunum$days), ]
jejunum_names <- rownames(jejunum)
aheatmap(c[jejunum_names, jejunum_names], Rowv = NA, Colv = NA, 
    cellwidth = 9, cellheight = 9, annRow = list(Day = as.factor(jejunum$days)), 
    breaks = 0.8)
```

<img src ="https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/jejunum%20correlation-1.png">

```r
# get summary statistic for every sample
kable(t(apply(c[jejunum_names, jejunum_names], 2, summary)), 
    format = "markdown")
```



|             |   Min.| 1st Qu.| Median|   Mean| 3rd Qu.| Max.|
|:------------|------:|-------:|------:|------:|-------:|----:|
|9980097039_B | 0.7889|  0.8905| 0.9182| 0.9090|  0.9361|    1|
|9980097056_C | 0.7625|  0.8717| 0.9044| 0.8980|  0.9254|    1|
|9982865068_D | 0.7586|  0.8747| 0.9103| 0.9022|  0.9327|    1|
|9982865068_E | 0.7774|  0.8897| 0.9212| 0.9116|  0.9394|    1|
|9372962074_G | 0.9029|  0.9358| 0.9641| 0.9538|  0.9703|    1|
|9372962074_L | 0.8569|  0.9375| 0.9541| 0.9486|  0.9687|    1|
|9377358063_B | 0.8611|  0.9306| 0.9550| 0.9485|  0.9682|    1|
|9377358066_A | 0.8883|  0.9358| 0.9580| 0.9488|  0.9684|    1|
|8909358002_D | 0.8778|  0.9335| 0.9591| 0.9519|  0.9712|    1|
|9287078038_K | 0.8449|  0.9380| 0.9512| 0.9477|  0.9703|    1|
|9287078041_D | 0.8654|  0.9192| 0.9328| 0.9279|  0.9430|    1|
|9287078048_J | 0.8493|  0.9357| 0.9500| 0.9462|  0.9684|    1|
|9287078051_L | 0.8907|  0.9321| 0.9623| 0.9527|  0.9706|    1|
|9377358019_L | 0.8853|  0.9425| 0.9627| 0.9524|  0.9689|    1|
|9377358053_D | 0.8862|  0.9408| 0.9624| 0.9558|  0.9742|    1|
|9377358065_J | 0.8961|  0.9260| 0.9611| 0.9513|  0.9700|    1|
|8909358005_F | 0.8632|  0.9303| 0.9455| 0.9429|  0.9665|    1|
|8909358006_F | 0.8264|  0.9099| 0.9382| 0.9273|  0.9556|    1|
|8909358007_A | 0.8621|  0.9367| 0.9548| 0.9434|  0.9632|    1|
|8909358007_L | 0.8736|  0.9338| 0.9556| 0.9446|  0.9645|    1|
|8909358275_C | 0.8697|  0.9384| 0.9496| 0.9424|  0.9588|    1|
|8909358279_K | 0.8593|  0.9335| 0.9566| 0.9430|  0.9619|    1|
|8981245016_C | 0.8041|  0.9040| 0.9150| 0.9124|  0.9307|    1|
|8981245022_J | 0.8903|  0.9318| 0.9541| 0.9456|  0.9614|    1|
|8981245037_B | 0.8603|  0.9377| 0.9490| 0.9441|  0.9618|    1|
|9031313009_A | 0.8785|  0.9370| 0.9634| 0.9540|  0.9699|    1|
|9031313023_A | 0.8972|  0.9369| 0.9578| 0.9490|  0.9624|    1|
|8909358007_C | 0.7858|  0.8790| 0.9155| 0.9068|  0.9394|    1|
|8909358279_G | 0.7935|  0.8851| 0.9189| 0.9089|  0.9377|    1|
|9031313012_H | 0.8818|  0.9341| 0.9541| 0.9462|  0.9655|    1|
|9031313015_H | 0.7586|  0.8586| 0.8907| 0.8791|  0.9105|    1|

No outliers to be removed as all samples have a mean correlation of at least 0.88 with the other jejunum samples.

### Blood and Lymph Nodes

```r
BLN <- c("Blood", "axillary_LN", "mesenteric_LN", "genital_pelvic_LN")
BLN <- MetaData[MetaData$tissue %in% BLN, ]
BLN <- droplevels(BLN)
BLN <- BLN[order(BLN$tissue, BLN$days), ]
BLN_names <- rownames(BLN)
aheatmap(c[BLN_names, BLN_names], Rowv = NA, Colv = NA, cellwidth = 2, 
    cellheight = 2, annRow = list(Tissue = BLN$tissue, Day = as.factor(BLN$days)), 
    breaks = 0.8)
```

<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/blood%20and%20lymph%20nodes%20correlation-1.png">

```r
# get summary statistic for every sample
kable(t(apply(c[BLN_names, BLN_names], 2, summary)), format = "markdown")
```



|             |   Min.| 1st Qu.| Median|   Mean| 3rd Qu.| Max.|
|:------------|------:|-------:|------:|------:|-------:|----:|
|9982865061_B | 0.8944|  0.9525| 0.9670| 0.9623|  0.9740|    1|
|9982865061_C | 0.9102|  0.9597| 0.9711| 0.9668|  0.9753|    1|
|9982865061_D | 0.9048|  0.9465| 0.9576| 0.9556|  0.9641|    1|
|9982865068_A | 0.9105|  0.9591| 0.9666| 0.9646|  0.9720|    1|
|9372962029_C | 0.9181|  0.9657| 0.9770| 0.9718|  0.9807|    1|
|9372962066_C | 0.9234|  0.9670| 0.9784| 0.9729|  0.9814|    1|
|9372962074_E | 0.9278|  0.9634| 0.9739| 0.9704|  0.9788|    1|
|9377361008_I | 0.9082|  0.9598| 0.9693| 0.9654|  0.9754|    1|
|8909358125_J | 0.9226|  0.9622| 0.9737| 0.9696|  0.9796|    1|
|9287078028_L | 0.9224|  0.9652| 0.9742| 0.9713|  0.9805|    1|
|9287078038_J | 0.9321|  0.9636| 0.9757| 0.9708|  0.9799|    1|
|9287078048_C | 0.9324|  0.9686| 0.9780| 0.9744|  0.9822|    1|
|9287078060_I | 0.9170|  0.9647| 0.9749| 0.9709|  0.9806|    1|
|9377358019_J | 0.9177|  0.9596| 0.9696| 0.9660|  0.9764|    1|
|9377358058_H | 0.9248|  0.9676| 0.9753| 0.9718|  0.9791|    1|
|9377358065_B | 0.9104|  0.9585| 0.9741| 0.9687|  0.9781|    1|
|9377361033_K | 0.9375|  0.9634| 0.9733| 0.9698|  0.9777|    1|
|8292044049_L | 0.9162|  0.9600| 0.9720| 0.9683|  0.9774|    1|
|8909358002_B | 0.9224|  0.9631| 0.9741| 0.9703|  0.9807|    1|
|8909358014_D | 0.9111|  0.9600| 0.9744| 0.9688|  0.9779|    1|
|8909358014_G | 0.9193|  0.9643| 0.9745| 0.9705|  0.9795|    1|
|8909358026_F | 0.9311|  0.9603| 0.9717| 0.9686|  0.9771|    1|
|8909358026_K | 0.9222|  0.9649| 0.9738| 0.9714|  0.9793|    1|
|8909358261_C | 0.9314|  0.9656| 0.9749| 0.9717|  0.9796|    1|
|8981245016_I | 0.8893|  0.9276| 0.9370| 0.9367|  0.9457|    1|
|8981245020_K | 0.9124|  0.9631| 0.9722| 0.9689|  0.9769|    1|
|8981245037_J | 0.9373|  0.9652| 0.9740| 0.9706|  0.9786|    1|
|9031313012_I | 0.9164|  0.9657| 0.9764| 0.9725|  0.9801|    1|
|9031313015_C | 0.9212|  0.9607| 0.9686| 0.9668|  0.9746|    1|
|9031313015_G | 0.9181|  0.9661| 0.9752| 0.9713|  0.9784|    1|
|8909358018_J | 0.9189|  0.9529| 0.9597| 0.9591|  0.9670|    1|
|8909358264_A | 0.9059|  0.9579| 0.9655| 0.9636|  0.9703|    1|
|8909358279_B | 0.9216|  0.9603| 0.9701| 0.9668|  0.9759|    1|
|9031313013_K | 0.9216|  0.9584| 0.9667| 0.9649|  0.9711|    1|
|9980097039_H | 0.8785|  0.9405| 0.9460| 0.9468|  0.9534|    1|
|9980097056_E | 0.8825|  0.9496| 0.9549| 0.9540|  0.9590|    1|
|9982865061_H | 0.8691|  0.9447| 0.9504| 0.9491|  0.9555|    1|
|9982865068_G | 0.8799|  0.9490| 0.9535| 0.9530|  0.9583|    1|
|9372962029_E | 0.9140|  0.9621| 0.9658| 0.9657|  0.9700|    1|
|9372962029_I | 0.9089|  0.9485| 0.9528| 0.9548|  0.9571|    1|
|9372962047_J | 0.9146|  0.9558| 0.9603| 0.9615|  0.9659|    1|
|9372962074_K | 0.8994|  0.9640| 0.9675| 0.9663|  0.9708|    1|
|9377358063_I | 0.8911|  0.9528| 0.9601| 0.9579|  0.9639|    1|
|8909358002_K | 0.9084|  0.9624| 0.9672| 0.9659|  0.9704|    1|
|9287078036_I | 0.9093|  0.9540| 0.9599| 0.9599|  0.9648|    1|
|9287078038_E | 0.8989|  0.9405| 0.9469| 0.9496|  0.9559|    1|
|9287078051_A | 0.9099|  0.9487| 0.9524| 0.9556|  0.9610|    1|
|9287078051_D | 0.9200|  0.9504| 0.9547| 0.9558|  0.9589|    1|
|9377358026_E | 0.9141|  0.9546| 0.9600| 0.9598|  0.9648|    1|
|9377361033_J | 0.9064|  0.9507| 0.9562| 0.9576|  0.9627|    1|
|9377361036_G | 0.9157|  0.9621| 0.9664| 0.9662|  0.9722|    1|
|8909358001_C | 0.9014|  0.9496| 0.9553| 0.9565|  0.9615|    1|
|8909358005_K | 0.8994|  0.9526| 0.9570| 0.9578|  0.9624|    1|
|8909358015_K | 0.8817|  0.9469| 0.9522| 0.9524|  0.9581|    1|
|8909358017_A | 0.9202|  0.9589| 0.9650| 0.9640|  0.9711|    1|
|8909358129_B | 0.9187|  0.9591| 0.9648| 0.9640|  0.9677|    1|
|8909358241_F | 0.9187|  0.9490| 0.9558| 0.9565|  0.9626|    1|
|8909358294_F | 0.9262|  0.9556| 0.9623| 0.9623|  0.9689|    1|
|8981245020_C | 0.9136|  0.9527| 0.9559| 0.9585|  0.9632|    1|
|8981245021_A | 0.9007|  0.9455| 0.9501| 0.9526|  0.9561|    1|
|8981245042_A | 0.9136|  0.9562| 0.9601| 0.9610|  0.9651|    1|
|9031313003_E | 0.9050|  0.9482| 0.9537| 0.9555|  0.9638|    1|
|9031313009_L | 0.8883|  0.9499| 0.9562| 0.9565|  0.9648|    1|
|8909358018_A | 0.9112|  0.9436| 0.9513| 0.9521|  0.9601|    1|
|8909358261_F | 0.8953|  0.9544| 0.9591| 0.9591|  0.9652|    1|
|8909358294_K | 0.8961|  0.9508| 0.9570| 0.9566|  0.9624|    1|
|9031313013_F | 0.8916|  0.9447| 0.9502| 0.9518|  0.9609|    1|
|9980097039_D | 0.9077|  0.9534| 0.9629| 0.9599|  0.9683|    1|
|9980097039_K | 0.9205|  0.9549| 0.9662| 0.9635|  0.9715|    1|
|9982865061_A | 0.8911|  0.9152| 0.9252| 0.9271|  0.9365|    1|
|9982865068_C | 0.9139|  0.9570| 0.9676| 0.9638|  0.9724|    1|
|9372962047_H | 0.9015|  0.9647| 0.9755| 0.9705|  0.9800|    1|
|9372962074_A | 0.9170|  0.9596| 0.9714| 0.9673|  0.9765|    1|
|9377361008_J | 0.9164|  0.9605| 0.9736| 0.9683|  0.9775|    1|
|9377361031_A | 0.9370|  0.9645| 0.9740| 0.9700|  0.9791|    1|
|8909358006_E | 0.9171|  0.9630| 0.9758| 0.9707|  0.9809|    1|
|9287078028_H | 0.9335|  0.9642| 0.9740| 0.9704|  0.9797|    1|
|9287078028_K | 0.9281|  0.9670| 0.9739| 0.9722|  0.9786|    1|
|9287078041_G | 0.9271|  0.9655| 0.9746| 0.9717|  0.9803|    1|
|9287078051_C | 0.9077|  0.9609| 0.9703| 0.9666|  0.9760|    1|
|9377358058_D | 0.9305|  0.9687| 0.9746| 0.9719|  0.9780|    1|
|9377361033_A | 0.9290|  0.9662| 0.9768| 0.9719|  0.9810|    1|
|9377361036_I | 0.9341|  0.9597| 0.9669| 0.9652|  0.9739|    1|
|9377361036_K | 0.9318|  0.9647| 0.9739| 0.9702|  0.9781|    1|
|8381688093_F | 0.9254|  0.9626| 0.9744| 0.9706|  0.9789|    1|
|8909358010_J | 0.9127|  0.9608| 0.9738| 0.9695|  0.9786|    1|
|8909358013_I | 0.9027|  0.9623| 0.9714| 0.9689|  0.9773|    1|
|8909358125_K | 0.9084|  0.9634| 0.9753| 0.9700|  0.9796|    1|
|8909358129_E | 0.9237|  0.9640| 0.9734| 0.9703|  0.9799|    1|
|8909358275_B | 0.9287|  0.9594| 0.9696| 0.9663|  0.9740|    1|
|8909358294_L | 0.9244|  0.9643| 0.9739| 0.9710|  0.9792|    1|
|8981245021_K | 0.9132|  0.9488| 0.9580| 0.9554|  0.9646|    1|
|8981245022_D | 0.9254|  0.9693| 0.9772| 0.9735|  0.9808|    1|
|8981245042_J | 0.9311|  0.9678| 0.9777| 0.9732|  0.9818|    1|
|9031313003_G | 0.9102|  0.9621| 0.9744| 0.9695|  0.9781|    1|
|9031313015_K | 0.9116|  0.9634| 0.9719| 0.9684|  0.9761|    1|
|8909358015_D | 0.9150|  0.9596| 0.9717| 0.9678|  0.9762|    1|
|8909358279_J | 0.9103|  0.9608| 0.9689| 0.9669|  0.9736|    1|
|8909358294_G | 0.9157|  0.9599| 0.9717| 0.9681|  0.9770|    1|
|9031313003_B | 0.9103|  0.9519| 0.9610| 0.9593|  0.9659|    1|
|9980097039_F | 0.9067|  0.9529| 0.9645| 0.9604|  0.9694|    1|
|9980097039_L | 0.9198|  0.9548| 0.9648| 0.9617|  0.9698|    1|
|9982865061_J | 0.9140|  0.9519| 0.9647| 0.9613|  0.9703|    1|
|9982865068_B | 0.8990|  0.9529| 0.9643| 0.9605|  0.9706|    1|
|9372962029_G | 0.9276|  0.9636| 0.9755| 0.9706|  0.9802|    1|
|9372962029_L | 0.9151|  0.9582| 0.9735| 0.9682|  0.9786|    1|
|9372962074_D | 0.9348|  0.9654| 0.9753| 0.9712|  0.9808|    1|
|9377358071_J | 0.9224|  0.9637| 0.9775| 0.9718|  0.9819|    1|
|8909358006_H | 0.9147|  0.9595| 0.9754| 0.9701|  0.9807|    1|
|9287078041_B | 0.9248|  0.9666| 0.9753| 0.9723|  0.9802|    1|
|9287078041_E | 0.9213|  0.9631| 0.9745| 0.9697|  0.9795|    1|
|9287078041_J | 0.9365|  0.9654| 0.9756| 0.9719|  0.9813|    1|
|9287078060_G | 0.9313|  0.9633| 0.9765| 0.9706|  0.9801|    1|
|9377358053_J | 0.9314|  0.9636| 0.9760| 0.9708|  0.9792|    1|
|9377358065_C | 0.9216|  0.9611| 0.9729| 0.9687|  0.9786|    1|
|9377361033_L | 0.9252|  0.9641| 0.9770| 0.9710|  0.9807|    1|
|8381688093_C | 0.9272|  0.9651| 0.9757| 0.9717|  0.9807|    1|
|8909358005_D | 0.9094|  0.9573| 0.9710| 0.9663|  0.9759|    1|
|8909358014_I | 0.9261|  0.9611| 0.9711| 0.9674|  0.9777|    1|
|8909358014_L | 0.9252|  0.9610| 0.9698| 0.9671|  0.9765|    1|
|8909358015_G | 0.9308|  0.9637| 0.9761| 0.9712|  0.9807|    1|
|8909358026_E | 0.9090|  0.9591| 0.9706| 0.9656|  0.9751|    1|
|8909358264_E | 0.9343|  0.9658| 0.9762| 0.9715|  0.9799|    1|
|8909358275_J | 0.9101|  0.9597| 0.9743| 0.9691|  0.9789|    1|
|8981245016_H | 0.8691|  0.9105| 0.9208| 0.9198|  0.9290|    1|
|8981245022_F | 0.9268|  0.9640| 0.9766| 0.9721|  0.9805|    1|
|8981245042_F | 0.9264|  0.9612| 0.9750| 0.9700|  0.9794|    1|
|9031313003_L | 0.9227|  0.9647| 0.9757| 0.9714|  0.9790|    1|
|9031313009_B | 0.9087|  0.9609| 0.9726| 0.9683|  0.9777|    1|
|8909358018_H | 0.9283|  0.9573| 0.9681| 0.9653|  0.9744|    1|
|8909358275_G | 0.9191|  0.9634| 0.9751| 0.9705|  0.9787|    1|
|8909358279_L | 0.9326|  0.9603| 0.9690| 0.9666|  0.9740|    1|
|9031313023_I | 0.9115|  0.9556| 0.9645| 0.9625|  0.9695|    1|

No outliers to be removed (mean correlation between blood/lymph samples of 0.879 or higher). We will therefore proceed to quantile normalization

## Quantile normalisation: 

```r
# perform quantile normalization
DATA <- normalize.quantiles(as.matrix(rawDATA), copy = FALSE)
DATA <- as.data.frame(DATA)
# log2 transform the data
DATA <- log2(DATA)

# Rename the rows and colunms
colnames(DATA) <- gsub("X", "", colnames(rawDATA))
row.names(DATA) <- row.names(rawDATA)

# display sample log2 intensity distribution
boxplot(DATA, range = 0, xaxt = "n", main = "sample log2 intensity distribution after quantile normalization")
```

<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/quantile%20normalization-1.png">

We can see, that after outlier removal and quantile normalization, the log2 intensity distribution looks good now.

## Export the new DATA and MetaDATA files:

```r
write.table(DATA, file = "DATA.txt", row.names = TRUE, col.names = NA)
write.table(MetaData, file = "MetaData_cleaned.txt", row.names = TRUE, 
    col.names = NA)
```
