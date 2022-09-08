# DCMI
## A Scalable Strategy for Accelerating Iterative Stencil Loops on FPGAs

### Table of contents
* [About](#About)
* [Technologies](#Technologies)
* [Setup](#setup)
* [Introduction](#Introduction)
* [Detailed Explanation](#Detailed-Explanation)
* [Single Time-Step](#Single-Time-Step)
* [Multiple Time-Stepsp](#Multiple-Time-Steps)

### About
This project introduces an architecture that optimizes Stencil Kernel Calculations and is offloaded to FPGAs with the use of Vitis (Vivado) High Level Synthesis tool (HLS).  

### Technologies
Project is created with:
* C/C+
* Vitis HLS 2021.2

### Setup
* Run Vitis (or Vivado) HLS adding the source & header files from the corresponding folder.
* The different no. of CKs is denoted as K.
* The no. of Time-Steps is denoted as D. 
* Each folder contains a version of the architecture with different factors K and D.
* The jacobi9d.cpp file should be used as the top function for the implementation of a single time-step.
* The SpatioTemporal.cpp file should be used as the top function for the implementation of multiple time-steps.
* Utilize the provided Test-Bench from the corresponding folder.
* In the header file the defined size of the grid can be modified.
* The number of Time-Steps can be modified in the Temporal.cpp file, by adding succesive calls to the jacobi9d funncion and declaring the intermediate variables.
* The number of CKs available requires modifications to the code of jacobi9d.cpp by adding/removing blocks that describe the Reuse Chains.

### Introduction 
