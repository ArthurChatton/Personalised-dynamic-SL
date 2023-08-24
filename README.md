# Personalised-dynamic-SL

R codes for the paper "Personalised dynamic super learning: an application in predicting convection volumes in hemodiafiltration" by Chatton et al.:
 - lrnr_XXX.r are functions adapted from the sl3 R package
 - SL_functions.r are the main functions implementing the method and the validation process
 - run_XXX.r are the master scripts related to the main analysis (XXX=continuous) and the sensitivity analysis (XXX=binary). The main analysis predicts a continous outcome, while the sensitivity analysis takes place in a classification setting.
 - visualisation.r is the script used to make the Figures and Tables.

For the setup used for the analysis, please see Supplementary Table 2 of the paper.
The seeds used are in the seeds.txt file.
