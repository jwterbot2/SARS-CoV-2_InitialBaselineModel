# SARS-CoV-2_InitialBaselineModel
Scripts used in [citation] to test parameter and model space for the construction of a baseline model of the intra-host 
population genetics of the SARS-CoV-2 virus. Reuse and modification of this code for use outside of this target organism
is encouraged with a request that this paper and/or repository are cited in any published material.

Briefly, the pipeline consists of:
0) Generating .ms files from simulations using the SLiM script.
1) Using the .ms files as input for the python script to calculate summary statistics.
2) Compiling, analyzing, and displaying the output from the python script using the R scripts.

More detailed readmes for each of these steps are included alongside the scripts in their respective folders.

Also included are example bash scripts which manage this pipeline for use with OSPool resources.
