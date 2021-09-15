Progress Report
================

Changes based on the final proposal
-----------------------------------
Our dataset and main objective have not changed since our proposal. Our main objective is to identify differentially expressed genes in seven different tissues in the first ten days after SIV infection. We have adapted the methods as described below to best answer our aims. The distribution of the work has remained the same. 

Progress of the analyses
------------------------

The data has been [imported, cleaned and normalized](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning.md): The tissues of interest were selected, reducing the sample size from 697 to 231. Three further samples were removed (outliers), resulting in a final number of 228 samples. These samples were quantile normalized and log2 transformed to allow statistical analysis.&nbsp;


**Question 1:** Is gene expression affected over time (i.e. by SIV infection) in different tissues?

A linear model for each of the seven tissues is fit to identify differential gene expression with respect to time in each tissue. As exploratory analysis we are testing a linear model with and without a quadratic term to see if a nonlinear model is a better fit to the data. 

1. Parameter specification: We were advised that including animal effect as a random effect and fitting a mixed effect model was out of the scope of this class, therefore when determining the parameters to include in the models we only had ‘time’ to consider. Therefore we selected ‘time’ and ‘time^2’ as the parameters for the linear models fit to each tissue. 

2. Fitting linear models: We fit 14 linear models (2 models per tissue type with one containing a quadratic of time and the other just looking at change in gene expression with respect to time) using the limma package. [Fitting linear models - tissue specific gene expression through time](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Exploratory%20Fitting/Exploration_of_different_models_for_tissue_specific_gene_expression.md)

3. Testing empirical reliability and goodness of fit: For each tissue, we will compare the linear model with the quadratic model, to determine which best represents the data. Using the F-tests in regression we will determine the Pbig (with quadratic) value and Psmall (without quadratic) value. By looking at the F distribution and  P-value we can plot a histogram and determine if the Pbig or Psmall value is preferred (i.e. if Pbig is significantly smaller than Psmall), meaning we will know if a quadratic function is preferred over the linear model for a given tissue type. 

4. Downstream Analysis: We have not yet proceeded to downstream analysis (GSEA), as we first need to determine which model fits the data better. 

**Question 2:** How do transcriptomic changes in blood compare to transcriptomic changes in other tissues?
Differentially expressed genes in blood that were identified by the linear model fit above will be compared to the genes identified in the other tissues. We will compare the genes by t-test rank.

**Question 3:** Are transcriptomic changes after SIV the same in the three different lymph nodes?

Thus far, 3 contrast matrices fitting linear models between lymph node type (genital-pelvic vs. axillary, genital pelvic vs. mesenteric, axillary vs. mesenteric) were generated. In terms of exploratory data analysis we will be accounting for the time and tissue effect on gene expression. For our current models, we have generated top genes that are differentially expressed, without time accounted for to just observe tissue effect between lymph nodes.
To analyze whether the transcriptomic changes after SIV infection are the same in the three different lymph nodes or whether they are different by location, we will perform unsupervised clustering for the lymph node samples to see whether the lymph nodes cluster by their location. If they do, we will generate a gene expression heatmap for the lymph node tissues containing all probes that are differentially expressed through time in at least one of the lymph nodes (based on the linear model developed in question 1. If the gene expression in the lymph node sites is similar, but temporally delayed (based on distance to infection site), we should see these patterns in the heatmap.
[Fitting linear models - lymph node contrast matrix](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%203/question3.md)

**Packages used so far in our analysis:** LIMMA (for the analysis of our microarray gene expression data), GEOquery (to help us get data from the NCBI Gene Expression Omnibus), Biobase, preprocessCore, NMF (to generate correlation heatmaps), Cluster, Ggplot2, and dplyr. 

1. Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016). Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. Annals of Applied Statistics 10(2), 946{963.
2. H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
 Davis, S. and Meltzer, P. S. GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics, 2007, 14, 1846-1847
3. Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M., Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M, MacDonald, J., Obenchain, V., Ole's, K. A, Pag'es, H., Reyes, A., Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan and M. (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” Nature Methods, 12(2), pp. 115–121. http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html. 
4. Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2016).  Cluster: Cluster Analysis Basics and Extensions. R package version 2.0.5.


Results 
-------

We have not been able to address all our aims yet, however we do have some primary results from fitting the linear models. When we fit a linear model for each tissue we see that there are differences in the number of differentially expressed genes between tissues (with an FDR cutoff of 0.05). This is also seen when we fit the linear model with the quadratic. When comparing the number of hits between the linear model with the quadratic to the linear model without the quadratic term, we see that including the quadratic results in a large reduction in the number of differentially expressed genes in the jejunum. Whether this reduction in the number of hits with the inclusion of the quadratic is due to increased fit of the model has yet to be determined. Interestingly, the inclusion of the quadratic term leads to more probes considered  to be differentially expressed in other tissues, e.g. tonsil. This suggests that the same model might not be best for all tissues. We need to test the goodness of fit of these models to determine which model fits the data best for each tissue before beginning downstream analysis. 
 
Based on our primary results, we observed that there are differentially expressed genes across different tissues at varying time points, which answers our first research question. We are considering the fact that preliminary analysis shows that there are differences in the differentially expressed genes between tissues as a positive result that supports our hypothesis. We are excited to continue to downstream analysis. 

&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Exploratory%20Fitting/Exploration_of_different_models_for_tissue_specific_gene_expression_files/figure-html/unnamed-chunk-5-1.png"> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Exploratory%20Fitting/Exploration_of_different_models_for_tissue_specific_gene_expression_files/figure-html/unnamed-chunk-5-2.png">


The above histograms show the p-value distribution for model 1 (without the quadratic) and model 2 (with the quadratic) of the Jejunum tissue. We can see from this distribution that model 1 may be a better fit for Jejunum due to the higher frequency of lower p-values. 

The number of genes differentially expressed for the tissue specific linear models can be seen in the table below: 

| Model Name| Parameters | Number of genes differentially expressed ( i.e. adjusted p-value <=0.05) |
| ------------- | ------------- | ------------- | 
| Model_Jejunum1 | days | 5598 | 
| Model_Jejunum2 | days+days^2 | 460 |  
| Model_Blood1 | days | 1756 | 
| Model_Blood2 | days+days^2 | 2827 |  
| Model_Tonsil1 | days | 100 | 
| Model_Tonsil2 | days+days^2 | 443 |  
| Model_ALN1 | days | 1858 | 
| Model_ALN2 | days+days^2 | 748 |  
| Model_MLN1 | days | 1666 | 
| Model_MLN2 | days+days^2 | 376 |  
| Model_GLN1 | days | 1245 | 
| Model_GLN2 | days+days^2 | 558 |  
| Model_Colon1 | days | 132 | 
| Model_Colon2 | days+days^2 | 386 |  



When fitting the linear models we ran into some issues when trying to determine whether or not to include the animal effect in the model. The original paper included the animal effect as a random effect and fit a mixed effect model. We came to the conclusion to exclude the animal effect because we are not including samples from the same animal in the same linear model (there’s only one sample/animal for each tissue type). Currently we are trying to determine the best way to assess the goodness of fit for the models. 
For fitting linear models for lymph nodes and generating contrast matrices, time was not taken into account for in the current model and so, it may be a better representation of differential expression of lymph nodes by remodelling to add the effect and interaction term between tissue*days.
The correlation analysis between transcriptomic changes in the blood and other tissues is not feasible with the methods taught in class. Therefore, we will initially simply compare the lists of differentially expressed genes between blood and other tissues by comparing the t-test ranking.
