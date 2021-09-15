Project Proposal
================

-   [Group information](#group-information)
-   [Introduction](#introduction)
-   [Division of labour](#division-of-labour)
-   [Dataset](#dataset)
-   [Aims and methodology](#aims-and-methodology)
-   [References](#references)

Group information
------------------

<table style="width:78%;">
<colgroup>
<col width="19%" />
<col width="19%" />
<col width="19%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th>Name</th>
<th>Department/Program</th>
<th>Expertise/Interests</th>
<th>GitHub ID</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Christina Michalski</td>
<td>Experimental Medicine (PhD)</td>
<td>Immunology</td>
<td><span class="citation">@ChristinaMi</span></td>
</tr>
<tr class="even">
<td>Heather A. MacKay</td>
<td>Wood Science (MSc)</td>
<td>Plant Biotechnology and molecular biology</td>
<td><span class="citation">@HAMacKay</span></td>
</tr>
<tr class="odd">
<td>Jeffrey Tang</td>
<td>Nondegree</td>
<td>Genomics and molecular biology</td>
<td><span class="citation">@jt1013</span></td>
</tr>
<tr class="even">
<td>Nichalle Brito</td>
<td>Biochemistry and molecular biology (MSc)</td>
<td>Proteomics and viral proteases</td>
<td><span class="citation">@Nichalle</span></td>
</tr>
<tr class="odd">
<td>Rebecca Asiimwe</td>
<td>Bioinformatics (MSc)</td>
<td>Computational Biology and Cancer Research</td>
<td><span class="citation">@rasiimwe</span></td>
</tr>
</tbody>
</table>

### Team Name: SIV in Rhesus Monkeys

Introduction
------------

Human immunodeficiency virus (HIV) has continued to be a public health issue on a global scale and since the start of the epidemic in 1981, about 70 million people have become infected and 35 million have died due to onset of acquired immunodeficiency syndrome (AIDS)-related diseases (1). Early diagnosis is critical as it enables patients to control for the infection via antiretroviral therapy in addition to reducing the risk of onward transmission (2). Despite our understanding that early detection can lead to prevention measures, it is estimated that half of the patients with HIV worldwide are “late-presenters” due to the asymptomatic nature of the virus (3). It is therefore evident that more knowledge is needed to understand the early stages of HIV.

Simian immunodeficiency virus (SIV), the source of HIV, does not affect humans but has a similar effect as HIV, in nonhuman primates (1,4). Thus, it is no surprise that SIV and HIV are closely related in viral replication and propagation. Rhesus monkeys susceptible to SIV are an ideal model for studying the early stages of HIV infection in humans due the similarity of SIV and HIV. In this project, we will analyze a publicly available dataset examining the initial stages (Day 0,1,3,7, and 10) of SIV infection (5). We hypothesize that temporal transcriptomic changes that are tissue-specific will be observed.

Division of labour
------------------

<table style="width:58%;">
<colgroup>
<col width="19%" />
<col width="19%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th>Project Job</th>
<th>Details</th>
<th>Name</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Literature Review</td>
<td>Pre and post data analysis literature review</td>
<td>All members</td>
</tr>
<tr class="even">
<td>Data cleaning</td>
<td> Import and data organization (1), background correction (2), and normalization (3)</td>
<td>Christina and Heather</td>
</tr>
<tr class="odd">
<td>Code writing</td>
<td> Parameter specification (1), and fitting the linear model (2)</td>
<td>Nichalle (1), Rebecca (2), Heather (1 and 2), Jeffrey (1,2)</td>
</tr>
<tr class="even">
<td>Validation</td>
<td>Testing empirical reliability and goodness of fit</td>
<td>Nichalle and Christina</td>
</tr>
<tr class="odd">
<td>Downstream Analysis</td>
<td>Gene enrichment analysis (1), Correlation analysis (2), and Principal Component Analysis (3) </td>
<td>Heather (1, 2) , Christina (2, 3), Rebecca (1, 3), Jeffrey (1, 2, 3)</td>
</tr>
<tr class="even">
<td>Deliverables</td>
<td>Poster and GitHub Repo</td>
<td>All members</td>
</tr>
</tbody>
</table>

Dataset
-------

Our dataset comes from a study examining the early host response to immunodeficiency virus infection using an SIV infection in rhesus monkey mode (5). Samples were obtained from 34 monkeys from 28 different tissues (8-25 samples per monkey), resulting in 697 samples that were analyzed by microarray (Illumina HumanHT-12 V4.0). On day 0, four monkeys were sacrificed (uninfected control) and the remaining animals were necropsied at day 1 (4 subjects), day 3 (9 subjects), day 7 (13 subjects) or day 10 (4 subjects). This allows analysis of tissue-specific and temporal gene expression changes in the first ten days of SIV infection. In the corresponding publication, tissue samples were grouped into four categories (blood, female reproductive tract, gastrointestinal tract and lymph nodes). We propose to perform more detailed analysis by distinguishing different tissue types. For more information about the data please see the [README.md](https://github.com/rasiimwe/Statistics-for-High-Dimensional-Biology---Project/tree/master/Data/Raw%20Data) file in the Raw Data folder of the repo. 

Aims and methodology
-------
### A) How does gene expression change in various tissues after SIV infection?

Analyze gene expression within the jejunum, blood, tonsil, axillary lymph node, mesenteric lymph node, genital pelvic lymph node, and colon on day 0 (control), 1, 3, 7, and 10. We excluded other tissue types due to the lack of at least 3 biological replicates for all time points.

#### Computational/statistical approaches:

Microarray data will be analyzed using RStudio statistical software and the Linear Models for Microarray Data (LIMMA) statistical package (Bioconductor version 3.4). To clean the data we will need to do a background correction and normalize the data. Then a linear model will be fit to the normalized data. Quality assessment and empirical reliability of the linear model will be tested to ensure accurate linear model fit. Once the design matrix has been specified the command &gt;eBayes will be used to determine the standard errors, a moderated t-statistic, and the log-fold-changes for each gene. To identify genes with tissue specific differential expression the &gt;decideTests, &gt;toptable and &gt;volcanoplot functions will be used. This methodology will be used to address all of the aims, however different parameters will be defined.

Specific for question A) Gene set enrichment analysis for all tissue types.

### B) Are there correlations between changes in blood transcriptomic profiles to tissues (i.e. identify markers in the blood that correlate to early responses in the tissues)?

There is a lot of variation between viral RNA on day 7 between different animals therefore it is possible to look at the correlation between blood markers and tissue specific gene expression changes.

#### Computational/statistical approaches:

A similar microarray analysis will be done for this question as stated above.

Specific for question B) Correlation analysis. Pearson correlation tests with a p-value less than or equal to 0.01 will be considered significant.

### C) Are there differences in the lymphatic response depending on the location of the lymph node (LN)?

We hypothesize that gene expression changes will be the same in all lymph nodes but be delayed based on their distance to the initial infection site (vagina). Therefore, we will analyze and compare the temporal transcriptomic changes between the genital-pelvic lymph-node (close to infection site), mesenteric lymph node (gut, medium distance from infection site) and the axillary lymph node (underarm, far from infection site).

#### Computational/statistical approaches:

A similar microarray and correlation analysis will be done for this question as stated above.

Specific for question C) PCA on the lymph node samples (i.e. We are predicting that time post infection to be a major principal component, but will the samples cluster by location in one of the major principal components?)

References
----------

1.  Stevens, D.R. et al. (2017). A Global Review of HIV Self-testing: Themes and Implications. AIDS Behavior. <doi:10.1007/s10461-017-1707-8>.
2.  Mocroft, A. et al. (2013). Risk factors and outcomes for late presentation for HIV-positive persons in Europe: results from the Collaboration of Observational HIV Epidemiological Research Europe Study (COHERE). PLoS Med 10: e1001510.
3.  Levy, I. et al (2016). Missed opportunities for earlier diagnosis of HIV in patients who presented with advanced HIV disease: a retrospective cohort study BMJ 6: e012721.
4.  Klatt, N.R., et al. (2012).Nonpathogenic Simian Immunodeficiency Virus Infections. Cold Spring Harbor Perspectives in Medicine 2: a007153.
5.  Barouch, D.H. et al. (2016). Rapid Inflammasome Activation Following Mucosal SIV infection of Rhesus Monkeys. Cell 165: 656-667.
