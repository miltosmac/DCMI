/***********************************
9-Point Jacobi 2-D  DCMI Header File
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

#define HEIGHT 160
#define WIDTH 160

#define D 3 //The number of time-steps that an acceleration scheme harvests parallelism from as the iteration depth.
#define R (WIDTH+1)//The maximum distance from a pattern element to the element at the center of the stencil.
#define I 1//Number of inputs to process concurrently
#define K (2*R*D+1) //Length of Reuse Buffer
#define N ((D + 4) * (D + 4)) //Length of Coeff Arrays
#define STREAM_DEPTH (WIDTH-(D+4)) //Depth of the FIFOs

#define MAX_ROWS 2160
#define MAX_COLUMNS 4096


typedef float mat;
//typedef ap_fixed<32,16> mat;

const extern mat coeffs_2D[81][N];

void Coefficient_Generation(
	);

void DCMI_Compute_Engine_Core(mat input,mat &output, int index);

void jacobi9d(
		hls::stream <mat> &in, //When you use an array as an argument to the top-level function. V. HLS assumes Memory is off-chip
		hls::stream <mat> &ut//Streams must be passed in and out of functions by-reference
		);
#endif // __JACOBI9TCAD__ not defined
