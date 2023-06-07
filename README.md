
# PowerCurvesASCA

This repository contains the code for reproducing the results of the paper Camacho J, Díaz C, Sánchez-Rovira P. Permutation tests for ASCA in multivariate longitudinal intervention studies. Journal of Chemometrics. 2022;e3398. doi:10.1002/cem.339816 

This repository requires the MEDA toolbox v1.3 at https://github.com/josecamachop/MEDA-Toolbox

Contact person: José Camacho (josecamacho@ugr.es)

Last modification of this document: 20/Jul/22


## Organization

The code is organized in a number of Matlab scripts:

- The script named "PowerCurves" perform the computation of complete power curves from a set of tests and plots results.

- The script named starting by "CreatePlot" plot the figures 

- The rest of scripts are named after the corresponding test in the paper. 

The present repository is organized as follows:

- The root folder contains anonymized data and the code for unconstrained permutations in the observations.
- html/ contains the publication of results with publish				

## How to use

To use this code with your own data, you need to substitute the data set that is imported using 'importdata' in the 'PowerCurves' files. E.g., at line 37 in PowerCurves_UncRAw.m. Please, start running the code from the root folder and then proceed with the factors and interaction.  