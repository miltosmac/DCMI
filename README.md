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

The architecture proposed in this section is greatly inspired by the work done in [[1]](#1). 
The proposed architecture of [[1]](#1) performs stencil calculation over an arbitrary number of 𝐷 time-domain
iterations. Essentially the design bypasses the calculation of intermediate time-steps and outputs results
directly for the target 𝐷.

Most of the ISLs’ core computation is a first degree, i.e., linear, polynomial. Constant coefficients are
multiplied with the data elements of the stencil pattern and the sum of these multiplications holds the
resulting value. These mathematical operations, i.e., addition and multiplication, are characterized by the
commutative and associative properties, as it follows, the computation of each term in the polynomial
can be calculated independently. Moreover, it is evident that the coefficients of the polynomial across
multiple iterations are only dependent on the stencil pattern and the depth of the iterations therefore
their calculation need not happen at runtime because the effective coefficients for every element can be
computed beforehand, in design time.

An architecture tailored to suite the calculation of Jacobi 9-Point Stencil is proposed. The design
utilizes the minimal amount of on-chip memory required to stream data and store them until all the
resulting data that are dependent on it have been calculated. Each element is fetched once from external
memory per accelerator invocation. The input and output of data transpires in a streaming manner and
lexicographical order is maintained.


### References
<a id="1">[1]</a> 
Mostafa Koraei, Omid Fatemi, and Magnus Jahre. 2019. DCMI: A Scalable Strategy for Accelerating Iterative Stencil Loops on FPGAs. ACM Trans. Archit. Code Optim. 16, 4, Article 36 (December 2019), 24 pages. https://doi.org/10.1145/3352813
