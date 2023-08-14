This part offers allows to calculate the optical constant of the ices of interest in ../lab_data/ and compare the with reference data from Pontoppidan et al. 2008.

By executing the run.sh-file in this directory, the optical constants (n/ k) of the laboratory ices are calculated and stored in ./results as well as for the follow-up step in ../5_optool/lnk_data

The code optical_const_calc_quick.py requests the user at each iteration to specify the parameters:
- Ice thickness (microns)
- Refractive index (no)
- Refractive index substrate (n2)
- MAPE ('mean absolute percentage error)

additionally, it is important to provide to following values for the 'optool'-calculations:
- Ice density (default: 0.47 um)
- Length of the wave number range in 1/cm (default: 1, 30)

This part of input parameters can be personalized/ unified for all ices from line 206 on.
Path and input specifications are possible there too.
