###############################################################################
README

TNER: Tri-Nucleotide Error Reducer

TNER is a novel background error suppression tool that provides a robust 
estimation of background noise to reduce sequencing errors using tri-nucleotide 
context data for better identification of low frequency somatic mutations in
ctDNA (Circulating tumor DNA).

Version: 1.1
Last updated: Sept/04/2018
Contact: shibing.deng {at} pfizer.com or tao.xie {at} pfizer.com / xietao2000 {at} gmail.com
Pfizer Early Clinical Development Biostatistics
Pfizer Early Oncology Development and Clinical Research
Copyright (c) 2018 Pfizer Inc.
###############################################################################

Pre-requisites/Installation:
----------------------------
Download and install R (version 3 or later)
Download and install PERL (version 5.8 or later)
Download and install the TNER packge and unzip the files

----------------------------
Test the TNER function using the demo dataset
----------------------------
"Rscript TNER_main.R TNER_example_input_file.txt hs_ave_bg_error.csv hs_depth.csv"

Note: “TNER_example_input_file.txt” is the example input data to be analyzed; the 2nd argument (default:
“hs_ave_bg_error.csv“) is a csv file with the average background error rate from healthy subjects; the
last argument (default: “hs_depth.csv”) is a csv file with the base coverage in healthy subjects.


----------------------------
Build background error profile using "Create_background_error_rate.R” 
----------------------------
"Rscript Create_background_error_rate.R"

Note: "Create_background_error_rate.R” creates an average background error rate file from individual healthy 
subject data (should be stored in the subfolder named as "healthy_subjects", file format is as same as 
"TNER_example_input_file.txt". The function outputs “hs_ave_bg_error.csv” and “hs_depth.csv” for the main 
function “TNER_main.R”


----------------------------
Process pileup files for TNER using pileup2actg.pl
----------------------------
"perl pileup2actg.pl 

Note: this perl script processes a pileup file to generate a tab-delimated txt file for function 
“TNER_main.R”.


----------------------------
Citation
----------------------------
TNER: A Novel Bayesian Background Error Suppression Method for Mutation Detection in Circulating Tumor DNA
Shibing Deng, Maruja Lira, Donghui Huang, Kai Wang, Crystal Valdez, Jennifer Kinong, Paul A. Rejto, Jadwiga R. Bienkowska, James Hardwick, Tao Xie
doi: https://doi.org/10.1101/214379 https://www.biorxiv.org/content/early/2017/11/05/214379
