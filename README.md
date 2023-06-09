
# PowerCurvesASCA

This repository contains the code for reproducing the results of the paper Camacho J, Díaz C, Sánchez-Rovira P. Permutation tests for ASCA in multivariate longitudinal intervention studies. Journal of Chemometrics. 2022;e3398. doi:10.1002/cem.339816 

This repository requires the MEDA toolbox v1.3 at https://github.com/josecamachop/MEDA-Toolbox

Contact person: José Camacho (josecamacho@ugr.es)

Last modification of this document: 20/Jul/22


## Organization

The code is organized in a number of Matlab scripts:

- The scripts named "PowerCurves" perform the computation of complete power curves from a set of tests and plots results.

- The scripts named starting by "CreatePlot" plot the figures in the paper (Note: to run "CreatePlotBoots.m" you need to download the aboxplot folder in the Github repo lemonzi/matlab at https://github.com/lemonzi/matlab)

- The rest of scripts are named after the corresponding test in the paper. 

The present repository is organized as follows:

- The root folder contains anonymized data and the code for unconstrained permutations in the observations.
- A/				Contains the code for tests for the treatment factor.
- AB/				Contains the code for tests for the interaction between treatment and time.
- B/				Contains the code for tests for the time factor.
- C(A)/				Contains the code for tests for the subject factor.
- Fig/				Contains the figures of the paper.
- html/ 				Contains the publication of selected results with publish.



## How to use

To use this code with your own data, you need to substitute the data set that is imported using 'importdata' in the 'PowerCurves' files. E.g., at line 37 in PowerCurves_UncRAw.m. Please, start running the code from the root folder and then proceed with the factors and interaction.  
