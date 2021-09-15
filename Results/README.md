# SIV in Rhesus Monkeys
---
## Results
---
## 1. Analyze how gene expression changes in various tissues after SIV infection

### A) Data Cleaning

Our dataset contained 697 samples from SIV infected Rhesus monkeys. After removal of tissues for which less than 3 replicates for each time point were available, 231 samples from 7 tissues remained (blood, colon, jejunum, tonsil as well as the axillary, mesenteric and genital-pelvic lymph node). For data cleaning, we downloaded the raw data files from the GEO website (GSE80013) and plotted correlation matrices. In general, correlation was high between all samples with higher correlation between samples from the same tissue (Fig.1A). To identify outliers, inter-sample correlation was analyzed within each tissue. 1 outlier from the colon (Fig. 1B) and two outliers from the tonsil were removed.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/correlation%20heatmap-1.png"> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
**Figure 1A: Inter-sample correlation before cleanup of all samples**

&nbsp;

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning_files/figure-html/colon%20correlation-1.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 1B: Inter-sample correlation before cleanup of colon samples**

### B) Fitting linear models 

For each tissue, we fit a linear model using time as a numeric variable. We compared linear models with and without a quadratic term for time. By plotting the top hits of the linear and the quadratic model, we can see that there are probes for which the quadratic or the linear model are a better fit (Fig.1C and Fig.1D).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Validation%20of%20Linear%20Models/Validation_of_Linear_Models_files/figure-html/unnamed-chunk-2-1.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 1C: Top 4 probes in the Jejunum identified by the simple linear model**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Validation%20of%20Linear%20Models/Validation_of_Linear_Models_files/figure-html/unnamed-chunk-2-2.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 1D: Top 4 probes in the Jejunum identified by the quadratic model**

To compare the two models in a more systematic manner, we looked at the difference in AIC for the linear and the quadratic model (Fig 1E and Table 1). For a very high percentage of probes, there is very little difference in AIC, indicating that there is no big advantage of choosing one model over the other. The "clear" cutoff at -2 is due to the way that the AIC is calculated: AIC = 2k - 2ln(L) where k is the number of model parameters and L is the maximized value of the likelihood function for the evaluated model. A difference of -2 simply reflects the higher number of degrees of freedom in the quadratic model when there is no difference in L for the two models.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Validation%20of%20Linear%20Models/Validation_of_Linear_Models_files/figure-html/unnamed-chunk-12-1.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 1E: Histogram of the difference in AIC values between the linear and the quadratic model for the jejunum**
&nbsp;
&nbsp;

**Table 1: Summary of 'AIC linear - AIC quadratic' for all probes in each tissue**

| | Jejunum |Blood | Tonsil | axillary_LN | mesenteric_LN | genital_pelvic_LN | Colon |
|--------|--------|--------|--------|--------|--------|--------|--------| 
|Min.|-2.00|-2.00|-2.00|-2.00|-2.00|-2.00|-2.00|
|1st Qu.|-1.82|-1.79|-1.80|-1.84|-1.84|-1.83|-1.78|
|Median|-1.23|-1.05|-1.12|-1.25|-1.29|-1.24|-1.02|
|Mean|-0.32|0.29|-0.14|-0.23|-0.39|-0.28|0.00|
|3rd Qu.|0.21|0.82|0.51|0.22|0.08|0.24|0.81|
|Max.|30.23|35.89|27.54|34.54|32.96|28.10|24.68|
|% of -2<=AIC.diff<=2|87.73|82.36|85.49|87.09|88.69|87.29|83.42|


There are no probes for which the linear model is significantly better but several probes for which the quadratic model is much better (AIC difference of up to 35). Considering moreover that the quadratic model includes a linear term, we will not lose any information when using the quadratic model. We therefore decide to use the quadratic model in all our (downstream) analyses.

The AIC comparison only tells other which of the two models is better. However, it doesn’t make any conclusion about the goodness of fit for each model. 

### C) Gene Enrichment Analysis

When applying the “quadratic” model, we identify several differentially expressed (DE) probes in all seven tissues (FDR <= 0.05) with significant differences in number of DE probes between tissues (Table 2).

**Table 2: Number of differentially expressed probes in each tissue**

|Tissue|Number of hits|
|--------|--------| 
|Jejunum|5463|
|Blood|4103|
|Tonsil|316|
|Axillary Lymph Node|2561|
|Mesenteric Lymph Node|2063|
|Genital Pelvic Lymph Node|1886|
|Colon|595|

We next asked which biological function these probes relate to and therefore conducted gene enrichment analysis of hits for each tissue via GO annotation (Table 3A), KEGG alignment (using “goana” and “kegga” from limma) as well as DAVID analysis. Interestingly, in all seven tissues, “ribosome” is the top associated KEGG term with the DE probes (Table 3B). Our results suggest that early changes in the host are mainly dominated by viral replication (for which the host ribosome’s are required), indicating that the virus is rapidly disseminating throughout the entire host.

**Table 3A: Annotated GO Terms for the Differentially Expressed Genes in the Mesenteric Lymph Node**

|GO Term|Description|
|--------|--------| 
|GO:0006613|cotranslational protein targeting to membrane|
|GO:0006614|SRP-dependent cotranslational protein targeting to membrane|
|GO:0045047|protein targeting to ER|
|GO:0072599|establishment of protein localization to endoplasmic reticulum|
|GO:0019058|viral life cycle|
|GO:0000184|nuclear-transcribed mRNA catabolic process, nonsense-mediated decay|
|GO:0009057|macromolecule catabolic process|
|GO:0006413|translational initiation|
|GO:0070972|protein localization to endoplasmic reticulum|
|GO:0006612|protein targeting to membrane|
|GO:1901566|organonitrogen compound biosynthetic process|
|GO:0006518|peptide metabolic process|
|GO:0043603|cellular amide metabolic process|
|GO:1901575|organic substance catabolic process|
|GO:0090150|establishment of protein localization to membrane|
|GO:0034655|nucleobase-containing compound catabolic process|
|GO:0072657|protein localization to membrane|
|GO:0019083|viral transcription|
|GO:0019439|aromatic compound catabolic process|
|GO:0046700|heterocycle catabolic process|


**Table 3B: Pathways that corresponds to the most enriched genes from the Mesenteric Lymph Node**

|Pathway|Identified Gene Count|
|--------|--------|
|Ribosome|28|
|Arrhythmogenic right ventricular cardiomyopathy (ARVC)|11|
|Cholinergic synapse|14|
|Ubiquitin mediated proteolysis|16|


## 2. Assess whether transcriptomic changes in blood are the same as in other tissues.
The number of DE probes range from 316 probes in the tonsil to 5463 probes identified in the jejunum. We asked whether there was some overlap between the DE probes from different tissues. As depicted in Figure 2A, there is considerable overlap between DE probes from the jejunum and the blood while there are also many probes that are DE only in one of the tissues. Similar results were obtained with all other tissues.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%202/question2_files/figure-html/unnamed-chunk-7-1.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 2A: Overlap of DE probes in blood and jejunum**

To more closely look at the “common” DE probes, we generated linear models comparing the time point against day 0 (--> four linear models: day1vs0, day3vs0, day7vs0, day10vs0) for each tissue separately and containing only probes that are DE in both blood and the tissue of interest (here: jejunum). The logFC between blood and jejunum were compared by fitting a linear regression curve (Fig. 2B and Table 4): There is correlation between the fold changes in the jejunum and in the blood, as evidenced by the R^2 values. The slope around 1 indicates that genes are up- and down-regulated to the same extent in blood and jejunum after SIV infection. At day 10 (slope = 1.4), the changes in jejunum are more pronounced than in blood for the genes that are DE in both tissues.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%202/question2_files/figure-html/unnamed-chunk-7-2.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 2B: Correlation between fold changes of common DE probes in the jejunum and blood**
&nbsp;
&nbsp;

**Table 4: Regression statistics from the correlation between fold changes of common DE probes in the jejunum and blood**

||Intercept| Slope | Adj. R^2|
|--------|--------|--------|--------| 
|Day 1 vs. 0|0.009| 1.10 | 0.57 |
|Day 3 vs. 0|0.0151| 0.92 | 0.55 |
|Day 7 vs. 0|0.0039| 1.11 | 0.63 |
|Day 10 vs. 0 |-0.0085| 1.40 |0.68 |

In conclusion, there is considerable overlap between DE expressed probes in a tissue with the DE probes in blood. These "common DE" probes are mainly affected in the same way in the tissues. However, one has to keep in mind, that there are many probes that are DE in a tissue-specific manner. Therefore, the transcriptomic changes upon SIV infection are tissue-specific while sharing common patterns across tissues. This approach is exploratory and does not compare the temporal gene expression profiles.

---
## 3. Investigate whether transcriptomic changes in lymph nodes vary by lymph node location.
We hypothesized that transcriptomic changes might be different in the three studied lymph nodes as distance from infection site varies. We find no evidence to support this hypothesis by using hierarchical clustering and PCA (Fig.3A and Fig.3b). 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%203/question3_files/figure-html/PCA-1.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 3A: PCA of the lymph node samples annotated by time point** 
&nbsp;
&nbsp;

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%203/question3_files/figure-html/PCA-2.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Figure 3B: PCA of the lymph node samples annotated by lymph tissue type**

The two plots are depicting the same clustering. In the first plot, the samples are annotated by their time point, in the second by their tissue type. We can see that the "day 0" samples cluster separately from the samples from the other days. The "day 0" samples are from uninfected control so it is not unsurprising that these samples are less similar to the other ones. Interestingly, there is no clustering by lymph node location, indicating that the transcriptomic changes between the lymph nodes are not significantly different.  In conclusion, we infer that the location of the lymph node and hence the distance from infection site does not have a significant impact on transcriptomic changes in the first 10 days after SIV infection. This suggests rapid viral dissemination throughout the entire host and/or rapid migration of immune cells that have encountered SIV antigens.

---
## Discussion and relation to published analysis 
A main limitation of this study is the study design: two thirds of the microarray samples had to be excluded from the analysis due to lack of replicates within that tissue/time group. Barouch et al. circumvented this problem by grouping samples into four tissue categories based on similarity (blood, lymphatic system, gastro-intestinal tract and female-reproductive tract). We find different probes to be DE in tissues from the same category (e.g. colon and jejunum). We did however not test whether these differences are statistically significant. In addition, based on our results from Gene Enrichment Analysis, although we identified plausible gene functions and found that the differentially expressed genes for these 7 tissue types all showed some relation to ribosomal proteins, it would be useful to provide further analyses to determine the directionality (low or high expression values) of expression for these genes involved in the suggested pathway to provide more insight toward making biological inferences.

Our analysis does not correct for viral load in each sample which was shown to differ between animals and between tissues within an animal by Barouch et al. and had an impact on (immune) gene expression.

It is also important to keep in mind that the macaque samples were analyzed on a human beadchip array. Despite >90% sequence similarity, transcriptomic differences are expected.

---
## Conclusion
In line with the findings from Barouch et al. we find that significant transcriptomic changes occur as early as day 1 post infection throughout host tissues. This knowledge is crucial for development of treatment as well as vaccines and more effective post-exposure prophylaxis drugs.



