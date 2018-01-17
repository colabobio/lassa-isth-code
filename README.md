## Lassa ISTH source code

This repository includes the R and Python scripts implementing the computational methods in the paper "Clinical and laboratory predictors of Lassa fever outcome in a dedicated treatment facility in Nigeria: an observational cohort study" by Peter Okokhere, Andres Colubri, et al.

There are four main Jupyter notebooks included, and one supporting notebook containing utilities:

1. Variable selection. This notebook implements the variable selection protocol. It starts with all variables in the data with P-value of univariate association with outcome lower than 0.1, and applied redundancy and cluster analyses using the [redun()](https://www.rdocumentation.org/packages/Hmisc/versions/4.1-0/topics/redun) and [varclus](https://www.rdocumentation.org/packages/Hmisc/versions/4.1-0/topics/varclus) functions in the Hmisc pacakge to remove superflous variables and eventually reach a parsimonious model.

2. Model fitting. This notebook uses multiple imputation with the [MICE package](https://cran.r-project.org/web/packages/mice/index.html), generating 100 multiple imputations for each model and fitting the model on each imputed copy of the dataset. The 100 parametrizations are then combined with [pool](https://www.rdocumentation.org/packages/mice/versions/2.46.0/topics/pool) function provided by MICE. Calculation of basic parameter measures of the models (AUC, adjusted R-squared, etc.) is done via bootstrap sampling with optimism correction.

3. Internal evaluation. The bootstrap data generated in the previous step is used in this notebook to carry out additional evaluation of the models, and to generate the plots that appear in Figure 2 in the paper.

4. Risk stratification. This notebook generates the plot that shows the mortality across each risk group (low, medium, high) as a function of the days of fever before presentation, which is a proxy for the delay in starting the treatment after symptom onset.

5. Logistic Regression Utils. This notebook is not meant to be run directly, but it is imported by the internal evaluation and risk stratification notebooks.
