#include "DCMI_Jacobi.h"

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)


void DCMI_Compute_Engine_Core(mat input,mat &output,int index){
#pragma HLS INLINE
	static ap_shift_reg<mat, STREAM_DEPTH> fifo_0,fifo_1,fifo_2,fifo_3,fifo_4,fifo_5;
	static mat RB[N]={0};//Result Buffer (storage for partial results)
#pragma HLS ARRAY_PARTITION variable=RB complete dim=1
	//for (int i=0;i<I;i++){//Load I data elements into buffer Input
		for(int j=0;j<N;j++){
			//RB[0] is only accessed when i=0 thus the result is not affected by the next (i++) iteration
			RB[j]+=coeffs_2D[index][j]*input;
	}
	output=RB[0];//Output Completed Entries
	Shift_RB_0:for (int c=0;c<(D+3);c++){//Shift the Result Buffer
		RB[c]=RB[c+1];
	}
	RB[D+3]=fifo_0.shift(RB[D+4],STREAM_DEPTH-1);
	Shift_RB_1:for (int c=(D+4);c<(2*D+7);c++){//Shift the Result Buffer - 0 to 2(D+4)-1
		RB[c]=RB[c+1];
	}
	RB[2*D+7]=fifo_1.shift(RB[2*D+8],STREAM_DEPTH-1);
	Shift_RB_2:for (int c=(2*D+8);c<(3*D+11);c++){//Shift the Result Buffer - 2(D+4) to 3(D+4)-1
		RB[c]=RB[c+1];
	}
	RB[3*D+11]=fifo_2.shift(RB[3*D+12],STREAM_DEPTH-1);
	Shift_RB_3:for (int c=(3*D+12);c<(4*D+15);c++){//Shift the Result Buffer - 3(D+4) to 4(D+4)-1
		RB[c]=RB[c+1];
	}
	RB[4*D+15]=fifo_3.shift(RB[4*D+16],STREAM_DEPTH-1);
	Shift_RB_4:for (int c=(4*D+16);c<(5*D+19);c++){//Shift the Result Buffer - 4(D+4) to 5(D+4)-1
		RB[c]=RB[c+1];
	}
	RB[5*D+19]=fifo_4.shift(RB[5*D+20],STREAM_DEPTH-1);
	Shift_RB_5:for (int c=(5*D+20);c<(6*D+23);c++){//Shift the Result Buffer - 5(D+4) to 6(D+4)-1
		RB[c]=RB[c+1];
	}
	RB[6*D+23]=fifo_5.shift(RB[6*D+24],STREAM_DEPTH-1);
	Shift_RB_6:for (int c=(6*D+24);c<(7*D+27);c++){//Shift the Result Buffer - 6(D+4) to 7(D+4)-1
		RB[c]=RB[c+1];
	}
	Add_Zero_Val:RB[N-1]=0;//Zero new values
}

void jacobi9d(
		hls::stream <mat> &in, //When you use an array as an argument to the top-level function. V. HLS assumes Memory is off-chip
		hls::stream <mat> &out//Streams must be passed in and out of functions by-reference
		)
{

#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=coeffs_2D
	int index;
	int counter=0;
	mat input;
	mat output;
	for (int i=0;i<HEIGHT+3;i++){
		for (int j=0;j<WIDTH;j+=I){
#pragma HLS PIPELINE II=1
				if (i < HEIGHT) {
					in.read_nb(input);
				}
				else { input = 0; }
				if(i==0 && j==0){
					index = 0;
				}
				else if(i==0 && j==1){
					index = 1;
				}
				else if(i==0 && j==2){
					index = 2;
				}
				else if (i == 0 && j == 3) {
					index = 3;
				}
				else if(i==0 && j>3 && j<WIDTH-4){
					index = 4;
				}
				else if (i == 0 && j == WIDTH - 4) {
					index = 5;
				}
				else if(i==0 && j==WIDTH-3){
					index = 6;
				}
				else if(i==0 && j==WIDTH-2){
					index = 7;
				}
				else if(i==0 && j==WIDTH-1){
					index = 8;
				}
				else if(i==1 && j==0){
					index = 9;
				}
				else if(i==1 && j==1){
					index = 10;
				}
				else if(i==1 && j==2){
					index = 11;
				}
				else if (i == 1 && j == 3) {
					index = 12;
				}
				else if(i==1 && j>3 && j<WIDTH-4){
					index = 13;
				}
				else if (i == 1 && j == WIDTH - 4) {
					index = 14;
				}
				else if(i==1 && j==WIDTH-3){
					index = 15;
				}
				else if(i==1 && j==WIDTH-2){
					index = 16;
				}
				else if(i==1 && j==WIDTH-1){
					index = 17;
				}
				else if(i==2 && j==0){
					index = 18;
				}
				else if(i==2 && j==1){
					index = 19;
				}
				else if(i==2 && j==2){
					index = 20;
				}
				else if (i == 2 && j == 3) {
					index = 21;
				}
				else if(i==2 && j>3 && j<WIDTH-4){
					index = 22;
				}
				else if (i == 2 && j == WIDTH - 4) {
					index = 23;
				}
				else if(i==2 && j==WIDTH-3){
					index = 24;
				}
				else if(i==2 && j==WIDTH-2){
					index = 25;
				}
				else if(i==2 && j==WIDTH-1){
					index = 26;
				}
				else if (i == 3 && j == 0) {
					index = 27;
				}
				else if (i == 3 && j == 1) {
					index = 28;
				}											  
				else if (i == 3 && j == 2) {						  
					index = 29;
				}													  
				else if (i == 3 && j == 3) {						  
					index = 30;
				}													  
				else if (i == 3 && j > 3 && j < WIDTH - 4) {		  
					index = 31;
				}													  
				else if (i == 3 && j == WIDTH - 4) {
					index = 32;
				}
				else if (i == 3 && j == WIDTH - 3){ 
					index = 33;
				}													  
				else if (i == 3 && j == WIDTH - 2) {
					index = 34;
				}
				else if (i == 3 && j == WIDTH - 1) {
					index = 35;
				}
				else if(i>3 && i<HEIGHT-4 && j==0){
					index = 36;
				}
				else if(i>3 && i<HEIGHT-4 && j==1){
					index=37;
				}
				else if(i>3 && i<HEIGHT-4 && j==2){
					index=38;
				}
				else if (i > 3 && i < HEIGHT - 4 && j == 3) {
					index = 39;
				}
				else if (i > 3 && i < HEIGHT - 4 && j == WIDTH - 4) {
					index = 41;
				}
				else if(i>3 && i<HEIGHT-4 && j==WIDTH-3){
					index=42;
				}
				else if(i>3 && i<HEIGHT-4 && j==WIDTH-2){
					index = 43;
				}
				else if(i>3 && i<HEIGHT-4 && j==WIDTH-1){
					index = 44;
				}
				else if(i==HEIGHT-1 && j==0){
					index = 72;
				}
				else if(i==HEIGHT-1 && j==1){
					index=73;
				}
				else if(i==HEIGHT-1 && j==2){
					index=74;
				}
				else if (i == HEIGHT - 1 && j == 3) {
					index = 75;
				}
				else if (i == HEIGHT - 2 && j == 0) {
					index = 63;
				}
				else if(i==HEIGHT-2 && j==1){
					index=64;
				}
				else if(i==HEIGHT-2 && j==2){
					index=65;
				}
				else if (i == HEIGHT - 2 && j == 3) {
					index=66;
				}
				else if(i==HEIGHT-3 && j==0){
					index=54;
				}
				else if(i==HEIGHT-3 && j==1){
					index=55;
				}
				else if(i==HEIGHT-3 && j==2){
					index=56;
				}
				else if (i == HEIGHT - 3 && j == 3) {
					index=57;
				}
				else if (i == HEIGHT - 4 && j == 0) {
					index=45;
				}
				else if (i == HEIGHT - 4 && j == 1) {
					index=46;
				}
				else if (i == HEIGHT - 4 && j == 2) {
					index=47;
				}
				else if (i == HEIGHT - 4 && j == 3) {
					index=48;
				}
				else if(i==HEIGHT-1 && j>3 && j<WIDTH-4){
					index=76;
				}
				else if(i==HEIGHT-2 && j>3 && j<WIDTH-4){
					index=67;
				}
				else if(i==HEIGHT-3 && j>3 && j<WIDTH-4){
					index=58;
				}
				else if (i == HEIGHT - 4 && j > 3 && j < WIDTH - 4) {
					index=49;
				}
				else if(i==HEIGHT-1 && j==WIDTH-1){
					index=80;
				}
				else if(i==HEIGHT-1 && j==WIDTH-2){
					index=79;
				}
				else if(i==HEIGHT-1 && j==WIDTH-3){
					index=78;
				}
				else if (i == HEIGHT - 1 && j == WIDTH - 4) {
					index=77;
				}
				else if(i==HEIGHT-2 && j==WIDTH-1){
					index=71;
				}
				else if(i==HEIGHT-2 && j==WIDTH-2){
					index=70;
				}
				else if(i==HEIGHT-2 && j==WIDTH-3){
					index = 69;
				}
				else if (i == HEIGHT - 2 && j == WIDTH - 4) {
					index = 68;
				}
				else if(i==HEIGHT-3 && j==WIDTH-1){
					index=62;
				}
				else if(i==HEIGHT-3 && j==WIDTH-2){
					index=61;
				}
				else if(i==HEIGHT-3 && j==WIDTH-3){
					index=60;
				}
				else if (i == HEIGHT - 3 && j == WIDTH - 4) {
					index=59;
				}
				else if (i == HEIGHT - 4 && j == WIDTH - 1) {
					index=53;
				}
				else if (i == HEIGHT - 4 && j == WIDTH - 2) {
					index=52;
				}
				else if (i == HEIGHT - 4 && j == WIDTH - 3) {
					index=51;
				}
				else if (i == HEIGHT - 4 && j == WIDTH - 4) {
					index=50;
				}
				else{
					index=40;
				}
				DCMI_Compute_Engine_Core(input, output, index);
				if (i * WIDTH + j > 3*WIDTH+2) { //D=2 so i-D and j+v-D would also be correct
					out.write_nb(output);
				}
		}
	}
	/*
	DCMI_Compute_Engine_Core( input, output, coeffs);
	out[HEIGHT - 1][WIDTH - 3] = output[0];
	DCMI_Compute_Engine_Core( input, output, coeffs);
	out[HEIGHT - 1][WIDTH - 2] = output[0];
	DCMI_Compute_Engine_Core( input, output, coeffs);
	out[HEIGHT - 1][WIDTH - 1] = output[0];
	*/
}
