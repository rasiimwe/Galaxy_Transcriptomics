# Evaluating Tissue- Specific and Temporal Gene Expression Changes with Simian Immunodeficiency Virus (SIV) in Rhesus Monkeys
---
## STAT 540 GROUP PROJECT
## TEAM: SIV-in-Rhesus-Monkeys
---
### Team Members:
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
<td>Experimental Medicine</td>
<td>Immunology</td>
<td><span class="citation">@ChristinaMi</span></td>
</tr>
<tr class="even">
<td>Heather A. MacKay</td>
<td>Wood Science</td>
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
<td>Biochemistry and molecular biology</td>
<td>Proteomics and viral proteases</td>
<td><span class="citation">@Nichalle</span></td>
</tr>
<tr class="odd">
<td>Rebecca Asiimwe</td>
<td>Bioinformatics</td>
<td>Computational Biology and Cancer Research</td>
<td><span class="citation">@rasiimwe</span></td>
</tr>
</tbody>
</table>

 ### Project Deliverables
- [Proposal](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Proposal/project_proposal.md)
- [Progress report](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Progress%20Report/progress_report.md)
- [Poster](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Poster/Poster-FINAL.pdf)
---
### How to load and run the code on your personal computer:
The code is written with the assumption that you clone the repository onto your desktop (i.e. Create project as subdirectory of: ~/Desktop) and the project directory is named "team_SIV-in-Rhesus-Monkeys". Additionally, because the data is too large to upload onto github, you MUST run the [Data_Cleaning.Rmd](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning.Rmd) to load the raw data and cleaned data to your personal device.

---
### 1.0 Introduction 
Human immunodeficiency virus (HIV) has continued to be a public health issue on a global scale and since the start of the epidemic in 1981, about 70 million people have become infected and 35 million have died due to onset of acquired immunodeficiency syndrome (AIDS)-related diseases [1](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/References.md). Early diagnosis is critical as it enables patients to control for the infection via antiretroviral therapy in addition to reducing the risk of onward transmission [2](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/References.md). Despite our understanding that early detection can lead to prevention measures, it is estimated that half of the patients with HIV worldwide are “late-presenters” due to the asymptomatic nature of the virus [3](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/References.md). It is therefore evident that more knowledge is needed to understand the early stages of HIV. Simian immunodeficiency virus (SIV), the source of HIV, does not affect humans but has a similar effect as HIV, in nonhuman primates [1,4](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/References.md). Thus, it is no surprise that SIV and HIV are closely related in viral replication and propagation. Rhesus monkeys susceptible to SIV are an ideal model for studying the early stages of HIV infection in humans due the similarity of SIV and HIV. To help us understand early host response to the virus, in this study we evaluated tissue-specific and temporal gene expression changes that occur after SIV infection in Rhesus Monkeys.

Our analysis was based on a dataset derived from a study that examined the early host response to immunodeficiency virus infection using an "SIV infection in rhesus monkeys" model.  This study was [published](https://www.ncbi.nlm.nih.gov/pubmed/27085913) in Cell on April 2016 and the corresponding dataset is publicly available on Gene Expression Omnibus under the accession number GSE80013. Figure 1 below depicts the study design that was used.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/Study%20Design.png">

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Figure 1: Study Design**

---
### 1.2 Hypothesis:
Tissue-specific transcriptomic changes occur in the first ten days after SIV infection in Rhesus Monkeys

---
### 1.3 Specific Objectives
1. Analyze how gene expression changes in various tissues after SIV infection
2. Assess whether transcriptomic changes in blood are the same as in other tissues.
3. Investigate whether transcriptomic changes in lymph nodes vary by lymph node location.

---
### Dataset used for analysis
Based on the dataset imported, 697 samples were analyzed (=697 columns) and for each each sample, 47,232 probes were quantified (=47,232 rows); a 697 x 47,232 matrix corresponding to over 30,000 genes.

---
### 2.0 Methods 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src = "https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/Methods.png"> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Figure 2: Methods applied for each of our objectives**

---
### 3.0 [Summary of Results](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/tree/master/Results)
### In detail:
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 1. [Data cleaning](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Data/Processed%20Data/Data_Cleaning.md) 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2. [Linear model fitting](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Exploratory%20Fitting/Exploration_of_different_models_for_tissue_specific_gene_expression.md)
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 3. [Validation of linear model fitting](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Validation%20of%20Linear%20Models/Validation_of_Linear_Models.md) 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 4. [Gene expression changes in various tissues after SIV infection](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/question1_pt1.md)
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 5. [Transcriptomic changes in blood the same as in other tissues](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%202/question2.md) 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 6. [Transcriptomic changes in lymph nodes vary by lymph node location](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%203/question3.md) 

---
### 4.0 Discussion
In our study we examined the early stages of host response to immunodeficiency virus infection using an SIV rhesus monkey model. The seven tissues for which the microarray data was analyzed, all show significant transcriptomic changes within the first 10 days after vaginal SIV infection. In general, the gene profiles are better described by a quadratic than a simple linear model which was validated by AIC comparison. Future analysis might include more sophisticated linear models (e.g. with time as factorial variable) to better describe the temporal pattern.

There are probes that are affected in a tissue-specific manner but also differentially expressed probes common to all tissues. These “common” probes tend to be affected in the same direction and to similar extent in blood and other tissues as determined by log fold change comparison. This approach is exploratory and does not compare the temporal gene expression profiles.
We hypothesized that transcriptomic changes might be different in the three studied lymph nodes as distance from infection site varies. We find no evidence to support this hypothesis by using hierarchical clustering, indicating a rapid viral dissemination and migration of immune cells in the lymphatic system. Other than Barouch et al., we did not group samples into four tissue categories. We find different probes to be DE in tissues from the same category (e.g. colon and jejunum). We did however not test whether these differences are statistically significant.

Interestingly, the main pathways affected relate to protein translation and transport, metabolic processes and viral life cycle in all studied tissues. This adds to the findings by Barouch et al. which focused on immune-related transcriptomic changes. Our results suggest that early changes in the host are mainly dominated by viral replication and dissemination. A more in-depth analysis is needed to confirm host immune response. The difference in results can be explained by different statistical approaches. Moreover, our analysis does not correct for viral load in each sample which was shown to differ between animals and between tissues within an animal by Barouch et al. and had an impact on (immune) gene expression.

A main limitation of this study is the study design: two thirds of the microarray samples had to be excluded from the analysis due to lack of replicates within that tissue/time group. It is also important to keep in mind that  the macaque samples were analyzed on a human beadchip array. Despite >90% sequence similarity, transcriptomic differences are expected 

---
### 5.0 Conclusion 
In line with the findings from Barouch et al. we find that significant transcriptomic changes occur as early as day 1 post infection throughout host tissues. This knowledge is crucial for development of treatment as well as vaccines and more effective post-exposure prophylaxis drugs.

---
### [References](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Report_Files/References.md)

