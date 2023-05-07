Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

# MATLAB code for analyzing tetrode, fiber photometry and behavioral data

## Description

This repository contains MATLAB code for analyzing data collected from mice performing a probabilistic Pavlovian conditioning task. The full analysis can be performed by running the main function called PV_pavlovian_analysis_main.m.

## Content

- .m files for data analysis
- license file

## Installation

- Download the .m file of this package and dependent packages (see dependencies below).
- Add directories (with subdirectories) with the downloaded .m files to the Matlab path.
- Download the sample data set from https://doi.org/10.6084/m9.figshare.22776776.v1.
- Run mountcb.m. Name the mounted Cellbase 'PV_pavlovian'. When prompted to 'locate CellBase database file', locate and choose the cellbase.mat file from the example data. Choose 'CellBase' for the time stamp conversion question.
- A typical installation takes about 10 minutes. 

## Demo

- Create a directory for the saved files.
- Download the example data set. We refer to the full path where the fiber photometry data is located on your computer as <data folder1>, and the path for the optogenetic inhibition data as <data folder2>.
- After installation, call PV_pavlovian_analysis_main(1, [1 1 1 1 1], [1 1], 1, [1 1 1 1 1 1 1 1 1], 1, [1 1], <data folder1>, <data folder2>) from the Command Window.
- For expected output, type help PV_pavlovian_analysis_main.
- Results will be saved in a 'data_analysis' folder in the CellBase directory of the sample data.
- Typical run time on demo data: 180 minutes.

## Instructions for use

- To run the data on your own CellBase, mount your Cellbase with mountcb.m
- In the 'Input argument check' panel of PV_pavlovian_analysis_main.m, replace 'PV_pavlovian' with the name of your CellBase
- Control which analysis panel is executed by the input arguments of PV_pavlovian_analysis_main (edit help PV_pavlovian_analysis_main for a full description on the input arguments)

## Dependencies

- This code relies on the CellBase data base system. MATLAB code for CellBase is available at www.github.com/hangyabalazs/CellBase.
- This code also depends on our MATLAB analysis package https://github.com/hangyabalazs/VP_data_analysis.
- Confidence intervals for regression are derived using the polypredci.m function (Star Strider, https://www.mathworks.com/matlabcentral/fileexchange/57630-polypredci, MATLAB Central File Exchange, retrieved December 30, 2020).

## System requirements  

Windows 10 64bits  
Any Intel or AMD x86-64 processor  
8 GB RAM  
MatlabR2016a or higher  
MATLAB Statistics Toolbox  

The code was also tested on Windows 10 64bits, Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz, 2208 MHz, 6 cores, 16 GB RAM, MATLAB R2020b.

Please contact us with any questions, bug reports, and general feedback:

**Panna Hegedus and Balazs Hangya**  
hangya.balazs@koki.hu
