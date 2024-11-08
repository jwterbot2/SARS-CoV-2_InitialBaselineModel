# SARS-CoV-2_InitialBaselineModel
Scripts used in Terbot et al 2023* to test parameter and model space for the construction of a baseline model of the intra-host 
population genetics of the SARS-CoV-2 virus. Reuse and modification of this code is encouraged with a request that this 
paper and/or repository are cited in any published material.

Briefly, the pipeline consists of:

0) Generating .ms files from simulations using the SLiM script.
1) Using the .ms files as input for the python script to calculate summary statistics.
2) Compiling, analyzing, and displaying the output from the python script using the R scripts.

More detailed readmes for each of these steps are included alongside the scripts in their respective folders.

Also included are example bash scripts which manage this pipeline for use with OSPool resources.

>*JW Terbot, BS Cooper, JM Good, JD Jensen. A Simulation Framework for Modeling the Within-Patient Evolutionary Dynamics of SARS-CoV-2. Genome Biology and Evolution. Volume 15:11. November 2023. https://doi.org/10.1093/gbe/evad204
