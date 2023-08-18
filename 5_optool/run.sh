#!/bin/bash

echo "Run this script to obtain optical constants and the opacity models of the ices via optool:"

echo "At first the optical constants (n-/ k-) of all ices in ./lab_data are determined and the results are written in the required .lnk-filkes"
bash ./4_optical_const/run.sh


#check on availabe data and do calculation for all of them
echo "Iterate through all available files in ../lab_data"
echo "With optical_const_calc_quick.py calculate the optical constantas (n/ k)"

echo "Run the script ..."
# run within the python environment
#python3 optical_const_calc_quick.py

echo "Compare to reference data and create plots"
# Scripts that return output in ./plots
#jupyter notebook nk_code_check.ipynb
