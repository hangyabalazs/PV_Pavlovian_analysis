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

Move the .m files on your MATLAB path. No further installation is needed.

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

Please contact us with any questions, bug reports, and general feedback:

**Panna Hegedus and Balazs Hangya**  
hangya.balazs@koki.hu
