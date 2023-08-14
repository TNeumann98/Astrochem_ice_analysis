#!/bin/bash

echo "Run this script to create all plots"

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist, create it!"
  mkdir plots
fi

#python3 from mymisc import *

#check on availabe data and do calculation for all of them
echo "Iterate through all available files in ../lab_data"
echo "With optical_const_calc_quick.py calculate the optical constantas (n/ k)"

echo "Run the script ..."
# run within the python environment
python3 optical_const_calc_quick.py

echo "Compare to reference data and create plots"
# Scripts that return output in ./plots
jupyter notebook nk_code_check.ipynb
