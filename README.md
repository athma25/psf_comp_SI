# Codes for solving plant-soil feedback with competition model
## Manusciprt titled "Coexistence of competing plants under plant-soil feedback" has the details of the model and results
## Authors: Athmanathan Senthilnathan, Rafael D'Andrea

## Files
* *parGen.sh*: Bash script to generate parameter files
* *nospace.m*: MATLAB function to numerically solve for a single parameter combination
* *d_eq.m*: MATLAB function for the ordinary differential equation system in Equation 1
* *nospace_sweep.sh*: Bash script to automate numerical solver for a parameter sweep
* *KmeansGapLib.R*: R library with functions for k-means clustering with gap statistic
* *nospace_summary.sh*: Bash script to automate k-means clustering for a parameter sweep
* *fig2.nb*: Mathematica notebook to produce figure 2
* *fig3.nb*: Mathematica notebook to produce figure 3, and Appendix A9 figures A9, A10, A11, A12, A13
* *A1.nb*: Mathematica notebook for Appendix A1
* *A5.nb*: Mathematica notebook for Appendix A5
* *A6.nb*: Mathematica notebook for Appendix A6
* *A7.nb*: Mathematica notebook for Appendix A7
* *A8.nb*: Mathematica notebook for Appendix A8 and Appendix A9 figure A14
* *myTxtFmt.m*: MATLAB script to change font size and interpreter
* *printPdf.m*: MATLAB script to print plot as PDF

## Instructions to perform a parameter sweep
1. Create *output* directory
2. Edit *parGen.sh* and run the bash script to create batch sub-directory under *output* and generate parameter files
> ./parGen.sh batch\_name
3. Run *nospace_sweep.sh* to numerically solve for each parameter combination for several random communities in the batch sub-directory. Code requires *batch_name*, *run_numbers* and *repeatition number*.
4. Run *nospace_summary.sh* to run k-means clustering for each repeatition of a parameter combination. Code requires *batch_name*, *run_numbers* and *repeatition number*.

## Codes tested in wsl2 with Ubuntu 20.04.6 and R3.6.3. MATLAB 2023a and Mathematica 13.3 installed in Windows 11
## Update: April 2024. Codes work in MATLAB 2022a and Mathematica 14.0
