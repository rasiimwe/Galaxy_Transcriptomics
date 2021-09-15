# Fitting tissue specific linear models and identifying DEG within each tissue: 
Group members responsible for this task: All members of team SIV in Rhesus Monkeys 

To fit linear models for each tissue we did the following: 
1. We conducted exploratory model fitting using limma and specifying different parameters (e.g. time as the only coefficient and time + time^2 as coefficients) 
    * For a more detained breakdown of exploratory model fitting see [Exploration_of_different_models_for_tissue_specific_gene_expression.md](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Exploratory%20Fitting/Exploration_of_different_models_for_tissue_specific_gene_expression.md)
    
2. Selected linear model based on AIC comparison
    * For a more detained breakdown of validation methodology see [Validation_of_Linear_Models.md](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models/Validation%20of%20Linear%20Models/Validation_of_Linear_Models.md)
    
To identify DEG within each tissue we did the following:
1. Generated a toptable for each tissue 
2. Filtered the toptable by FDR <= 0.05
    
See [question1_pt1.md](https://github.com/STAT540-UBC/team_SIV-in-Rhesus-Monkeys/blob/master/Methods/Question%201/Linear%20Models%20/question1_pt1.md) for full methodology and outputs.
