/***********************************
9-Point Jacobi 2-D  DCNI Header File
***********************************/
#include <stdio.h>
#include "hls_stream.h"
#include "ap_shift_reg.h"
#include <assert.h>
#include <ap_fixed.h>

// Compare TB vs HW C-model and/or RTL
#ifndef __JACOBI9TCAD__
#define __JACOBI9TCAD__

#define HW_COSIM


#define HEIGHT 2560
#define WIDTH 2560

#define D 2 //The number of time-steps that an acceleration scheme harvests parallelism from as the iteration depth.
#define R (WIDTH+1)//The maximum distance from a pattern element to the element at the center of the stencil.
#define I 1//Number of inputs to process concurrently
#define K (2*R*D+1) //Length of Reuse Buffer
#define N ((D + 3) * (D + 3)) //Length of Coeff Arrays
#define E (2*R*D+I) //Length of Reuse Buffer in case of Spatial Parallelism
#define STREAM_DEPTH (WIDTH-(D+3)) //Depth of the FIFOs

#define MAX_ROWS 2160
#define MAX_COLUNNS 4096

//typedef float mat;
typedef ap_fixed<32,16> mat;

//GLOBAL VARIABLES
/*
extern mat coef0_0[N], coef0_1[N], coef0_2[N], coef1_0[N], coef1_1[N], coef1_2[N], coef2_0[N], coef2_1[N], coef2_2[N], coef0_3W4[N], coef1_3W4[N], coef2_3W4[N];
extern mat	coef0_W1[N], coef0_W2[N], coef0_W3[N], coef1_W1[N], coef1_W2[N], coef1_W3[N], coef2_W1[N], coef2_W2[N], coef2_W3[N], coef3H4_W1[N], coef3H4_W2[N], coef3H4_W3[N];
extern mat	coefH1_0[N], coefH1_1[N], coefH1_2[N], coefH2_0[N], coefH2_1[N], coefH2_2[N], coefH3_0[N], coefH3_1[N], coefH3_2[N], coef3H4_0[N], coef3H4_1[N], coef3H4_2[N];
extern mat	coefH1_W1[N], coefH1_W2[N], coefH1_W3[N], coefH2_W1[N], coefH2_W2[N], coefH2_W3[N], coefH3_W1[N], coefH3_W2[N], coefH3_W3[N], coefH3_3W4[N], coefH2_3W4[N], coefH1_3W4[N], coeffs[N];
*/
const extern mat coeffs_2D[49][N];


void mod_div(int dividend, int& remainder, int& quotient);

void Coefficient_Generation();

void DCMI_Compute_Engine_Core(mat input,mat &output,int index);

void jacobi9d(
		hls::stream <mat> &in, //When you use an array as an argument to the top-level function. V. HLS assumes Memory is off-chip
		hls::stream <mat> &out //Streams must be passed in and out of functions by-reference
		);
#endif // __JACOBI9TCAD__ not defined
