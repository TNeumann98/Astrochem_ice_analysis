# Astrochem_ice_analysis
This repository should provide an extended tool for the analysis of spectral data and the investigation of ices (frozen molecules) in proto-stellar envelopes and other cool environments. 

The aim is the unification and optimization of pre-existing codes to create a unified routine to make radiative transfer modelling possible.

Input files: [add manually to the respective directory: /lab_data, /obs_data, /ref_data]
- laboratory absorbance spectra of molecules at specific temperatures
  - data can be found in the LIDA-database [https://icedb.strw.leidenuniv.nl/]
  -  the folder structure is adjusted to the format in which files can be downloaded from LIDA
- flux of observational spectrum of the desired object
- for reference: the absorbance spectra in the same form as the lab spectra
  
Output in different steps:
- optical constants (n/k) for laboratory data
- optical depth for molecules of choice
- RADMC-3D modelling of selected source with the chosen laboratory data

Prequisitions:
Included codes which require separate activation:
- optool
- RADMC-3D

the procedure for these codes can be found in the respective README-files in the folders 5/6.

General user manual:
1.) $git clone$ this repository to your local machine
2.) Activate the codes mentioned above
3.) Search and add the input files you wish to investigate
4.) Adjust the variables in the top-most parameter.py file to your setup of interest
5.) Open a terminal and run
$python3 run.sh$ : this file executes the scripts (run.sh) in each folder which run the respective .py-files which are also creating the plots which can be found in the plots-folder. 
