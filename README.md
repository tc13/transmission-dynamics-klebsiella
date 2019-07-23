# transmission-dynamics-klebsiella
Scripts to allow replication of the analysis from "Transmission dynamics and control of multidrug-resistant Klebsiella pneumoniae in neonates in a developing country"

The "risk factor models" can be run by following the risk_factor_models_klebsiella.R script in R-studio. This script will load in data for each patient day at risk of becoming colonised from Data/klebsiella_acquisition_reproducible. The Bayesian statistical models are run by calling scripts from the Models_Stan folder. 

The "transmission models" which incorporate genomic data in the form of sequence typying isolates of ESBL Klebsiella pneumoniae s.l.. These are run by following the code in transmission_models_klebsiella.R in R-Studio.

It is important that all required packages are installed. There are notes on this at the top of each script. In particular please ensure that rstan is correctly installed and configured for your system. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details. 

For any comments on this code, please contact me on thomas.crellen@ndm.ox.ac.uk or tomcrellen@gmail.com. The code is my own, the original dataset is the property of Prof Ben Cooper, Prof Paul Turner and Dr Claudia Turner.
