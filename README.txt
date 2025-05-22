# SRIR Extrapolation via Image–Source Modelling & Sub-space Decomposition  
*Code accompanying* **“Spatial room impulse response extrapolation using the image-source method”** (LAC-25, 2025)

---

## Overview
This repository contains a MATLAB reference implementation of the method proposed in the paper above.  
Starting from **one** higher-order Ambisonic spatial room impulse response (SRIR) and a shoebox geometry, the algorithm

1. decomposes the SRIR into **salient** (direct + early reflections) and **diffuse** sub-spaces using GSVD,  
2. predicts the direct sound and early reflections at arbitrary source/receiver positions with the **image-source model** (ISM),  
3. maps each salient arrival to the target position by *delay*, *gain* and *HOA rotation*, and  
4. re-assembles the extrapolated salient part with the unaltered diffuse tail.

The code reproduces all quantitative results and figures in the paper.

---


Download all files.
Run Matlab Testing code in \FinalCode\UseSingleRecordings.m .

* now can download this file automaticly with matlab code
Download the file "6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa 564.6 MB"  from https://zenodo.org/records/5720724.
Put it in to the fold: " 6dof_SRIRs_eigenmike_SH "


Listening Test Audio File in \listenTest.
Except coding in \FinalCode, other coding and file is from opensouce toolbox. put in this just for can easily run test code.
