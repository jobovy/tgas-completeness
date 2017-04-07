# tgas-completeness

The Gaia DR1 TGAS completeness and a new stellar inventory of the solar neighborhood

## Overview

This code repository contains all of the code to reproduce the results
from the paper of Bovy (2017, in prep.) on determining the Gaia DR1
TGAS selection function and using it to determine TGAS' completeness
for different stellar types and for measuring the local number density
and vertical stellar density profile of different types of
main-sequence and giant stars.

## AUTHOR

Jo Bovy - bovy at astro dot utoronto dot ca


## Code

## 1. [py/TGAS-selection-function.ipynb](py/TGAS-selection-function.ipynb)

(render this notebook on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/tgas-completeness/blob/master/py/TGAS-selection-function.ipynb), where you can toggle the code)

This notebook contains all of the code to determine the raw TGAS
selection function as a function of (*J*,*J-Ks*,RA,Dec), as discussed
in Appendix A and B.

## 2. [py/TGAS-effective-completeness.ipynb](py/TGAS-effective-completeness.ipynb)

(render this notebook on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/tgas-completeness/blob/master/py/TGAS-effective-completeness.ipynb), where you can toggle the code)

This notebook contains the various calculations of the effective
completeness and the effective volume completeness described in
Section 3 of the paper.

## 3. [py/TGAS-stellar-densities.ipynb](py/TGAS-stellar-densities.ipynb)

(render this notebook on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/tgas-completeness/blob/master/py/TGAS-stellar-densities.ipynb), where you can toggle the code)

This notebook holds all of the code that determines the stellar
density profiles of different stellar subtypes along the main sequence
and along the giant branch. All of the analysis in Sections 4, 5, and
6 of the paper is contained in this notebook.

## 4. [py/effsel.py](py/effsel.py)

This badly named python module contains the definitions of the main
sequence and giant branch as used in the paper.
