We moved to https://github.com/SayakaMiura/CloneFinder.

CloneFinder_v0.1.1
Updated June 6, 2020
==================

CloneFinder was developed by Sudhir Kumar

CloneFinder is a method aimed at estimating individual clone genotypes and frequencies within a tumor sample using a phylogenetic approach.  

Installation
==================
CloneFinder is a python script developed in a Windows and Unix/Linux 64-bit architecture, and does not need to be compiled prior to usage. To use CloneFinder, the following python packages need to be installed.

Python 3
SciPy
NumPy
Biopython 

Additionally, CloneFinder uses MEGA-CC (>= 7.0.18). MEGA-CC can be downloaded for Windows, Mac OS X, and Linux from (http://www.megasoftware.net). SciPy, NumPy, and Biopython are python packages and can be installed by using Anaconda (https://www.anaconda.com/) or PyPI (https://pypi.org/).
CloneFinder for python 2 is also available. However, this version has not been maintained, and bugs are not fixed.  

Input file
==================
The input file is a tab-delimited text file, which contains the read counts of wild-type and mutant alleles for each sample. Normal sample should not be included. Each line in the input files gives information for each variant. Explanation of each column and specific formatting can be found below.
 
* "XX:ref": Reference read count for the sample, XX
* "XX:alt": Variant read count for the sample, XX


***IMPORTANT***
1. Low quality SNVs (e.g., low coverage) should be removed from the Input file.
2. The data needs to contain parsimony-informative sites, when the tumor profiles of presence/absence of mutations are generated. The recommended number of samples per dataset is more than 4. 
3. Input file name and tumor/sector names can contain only alphabetic and numeric characters and should not begin with numbers. Total read count for each variant cannot contain values  < 1. 
4. cancer cell fraction (CCF) data needs to be converted into SNV read count table, under the assumption of CNA-free SNVs. Thus, variant read count is CCF that is multiplied with the total read count and divided by two. Reference read count can be computed by subtracting variant read count from total read count. 


Parameters
==================
The following parameters can be changed by editing CloneFinder/options.ini file.
* FreqCutoff
Cutoff of estimated clone frequencies. When estimated clone frequencies are less than this Cutoff value, estimated clone frequencies are reported as zero (Default: 0.02).
* Total_Read_Count_CutOff
Cutoff of total read count. Clone frequencies are estimated by using SNVs that have larger number of this cutoff(Default: 50). 
* Mutant_Read_Count_CutOff
Cutoff of variant read count. Clone frequencies are estimated by using SNVs that have larger number of this cutoff (Default: 2). 


Example 
==================
An example dataset (Example_data/Input.txt) is provided to run CloneFinder with default parameters (see below for individual parameter defaults). To run CloneFinder on InputTest.txt, please follow commands below from the CloneFinder directory.

       python clonefinder.py snv path/Input.txt

After running CloneFinder, the following output files can be found in the directory of input file. 


Output files
==================
1. Clone genotypes (e.g., Input_CloneFinder.meg)
This file can be opened with a text editor or MEGA (http://www.megasoftware.net). In each clone sequence, 'A' indicates wild-type allele and 'T' indicates mutant-type allele. The order follows that of the original input file. 
2. Clone frequencies (e.g., Input_CloneFinder.txt)
Clone frequencies in each tumor sample are given. 
3. Clone phylogeny (e.g., Input_CloneFinder.nwk)
This file can be opened with a text editor or MEGA (http://www.megasoftware.net). 
4. Summary table (e.g., Inputsnv_summary.txt)
Parameter settings used for the analysis are listed. 

How to cite
=================
If you use this CloneFinder software in your work, please cite the accompanying publication:

Sayaka Miura, Karen Gomez, Oscar Murillo, Louise A Huuki, Tracy Vu, Tiffany Buturla, & Sudhir Kumar. Predicting clone genotypes from tumor bulk sequencing of multiple samples. Bioinformatics (2018) 34(23):4017-4026
