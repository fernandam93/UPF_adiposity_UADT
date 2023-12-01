# UPF_adiposity_UADT

This code relates to the project titled "Ultra‐processed foods, adiposity and risk of head and neck cancer and oesophageal adenocarcinoma in the European Prospective Investigation into Cancer and Nutrition study: a mediation analysis" published in the European Journal of Nutrition (https://doi.org/10.1007/s00394-023-03270-1).

## Data management:

**First, we created a new dataset including the variables that would be required in all the analyses to avoid having to create the same variables multiple times.**

00-data_management.R

## Descriptive characteristics:

01-descriptive.R

## Main and secondary analyses:

02-Cox_regressions_N4.R         
> **_This file also includes sensitivity analyses in g/d_**

## Mediation analysis:

03-exploring_mediators_N4_cont.R

04-mediation_w_interactions_N4_bootstrapping_directCF_cont.R

## Sensitivity analyses:

05-Cox_regressions_N4_water_adjusted.R

06-Cox_regressions_N4_energy_adjusted.R

07-Cox_regressions_N4_exc_all_2yrs.R

08-Cox_regressions_N4_energy.R

09-Cox_regressions_N4_cc_analysis.R

10-Cox_regressions_N4_multiple_imputation.R

11-Cox_regressions_Nall_deaths.R

## Forest plots:

12-main_forest_plots_onlyN4.R

13-subtypes_forest_plots_onlyN4.R

14-main_forest_plots_sensitivity_kcal_g_relANDabs_onlyN4.R

## Subtype heterogeneity analysis:

15-subtypes_heterogeneity.R

## Interactions and stratified analyses:

16-Interactions_alc_AC.R

17-Interactions_alc_HNC.R

18-Interactions_education_AC.R

19-Interactions_education_HNC.R

20-Interactions_PA_AC.R

21-Interactions_PA_HNC.R

22-Interactions_sex_AC.R

23-Interactions_sex_HNC.R

24-Interactions_smoking_AC.R

25-Interactions_smoking_HNC.R

## Assumptions:

26-assumptions.R
