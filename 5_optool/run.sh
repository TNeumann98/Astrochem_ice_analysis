#!/bin/bash

echo "Run this script to obtain the opacity models of the ices via optool:"

echo "At first the optical constants (n-/ k-) of all ices in ./lab_data are determined and the results are written as the required .lnk-files in ./lnk_data"
echo "Check if optool is properly installed!"

#check on availabe data and do calculation for all of them
echo "Iterate through all available files in ./lnk_data"

if [ ! -d "plots" ]; then
  echo "Directory does not exist, create it!"
  mkdir plots
fi
echo "Compare to reference data and create plots"
# Scripts that return output in ./plots
#jupyter notebook Plot_opt_const.ipynb

#specify the ices
ice=( silicate pureCO210K 100H2O10CO210K 100CO70CO210K2 100CO26CO210K2 100CO4CO210K2 )
#files to run
nice="d0.00025_n01.16_n21.73_err0.0001_pureCO2_10K_test_quick5_kn_constants"
echo "Iterate through the available files, calculating the following:"

#for ice in "${ices[@]}";do
#    echo "$ice"
if [ ! -d "results/" ]; then
  mkdir "results"
fi


echo "Run the script ..."
./optool pyr-mg70 0.87 c 0.13 -l 0.1 1e4 10000 -a .1 1. -p 0.25 -m lnk_data/"$nice".lnk 0.7 -cde -radmc
cp dustkappa.inp ./"results/dustkappa_ "$nice".inp"
