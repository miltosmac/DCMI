#include "DCMI_Jacobi.h"

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)


void DCMI_Compute_Engine_Core(mat input,mat &output,int index){
#pragma HLS INLINE
	//static hls::stream<mat> fifo_0 ("RB_FIFO_0"),fifo_1 ("RB_FIFO_1"),fifo_2 ("RB_FIFO_2"),fifo_3 ("RB_FIFO_3");
	static ap_shift_reg<mat, STREAM_DEPTH> fifo_0,fifo_1,fifo_2,fifo_3;
	static mat RB[N]={0};//Result Buffer (storage for partial results)
/*DO_PRAGMA(HLS stream depth=STREAM_DEPTH variable=fifo_0)
DO_PRAGMA(HLS stream depth=STREAM_DEPTH variable=fifo_1)
DO_PRAGMA(HLS stream depth=STREAM_DEPTH variable=fifo_2)
DO_PRAGMA(HLS stream depth=STREAM_DEPTH variable=fifo_3)*/
#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=RB
	Mult_Acc_Circuit:for(int j=0;j<N;j++){
		RB[j]+=coeffs_2D[index][j]*input;
	}
	output=RB[0];//Output Completed Entries
	Shift_RB_0:for (int c=0;c<(D+2);c++){//Shift the Result Buffer
		RB[c]=RB[c+1];
	}
	RB[D+2]=fifo_0.shift(RB[D+3],STREAM_DEPTH-1);
	Shift_RB_1:for (int c=(D+3);c<(2*D+5);c++){//Shift the Result Buffer - (D+3) to 2(D+3)-1
			RB[c]=RB[c+1];
	}
	RB[2*D+5]=fifo_1.shift(RB[2*D+6],STREAM_DEPTH-1);
	Shift_RB_2:for (int c=(2*D+6);c<(3*D+8);c++){//Shift the Result Buffer - 2(D+3) to 3(D+3)-1
			RB[c]=RB[c+1];
	}
	RB[3*D+8]=fifo_2.shift(RB[3*D+9],STREAM_DEPTH-1);
	Shift_RB_3:for (int c=(3*D+9);c<(4*D+11);c++){//Shift the Result Buffer - 3(D+3) to 4(D+3)-1
			RB[c]=RB[c+1];
	}
	RB[4*D+11]=fifo_3.shift(RB[4*D+12],STREAM_DEPTH-1);
	Shift_RB_4:for (int c=(4*D+12);c<(5*D+14);c++){//Shift the Result Buffer - 4(D+3) to 5(D+3)-1
			RB[c]=RB[c+1];
	}
	Add_Zero_Val:RB[N-1]=0;//Zero new values
}

void jacobi9d(
		hls::stream <mat> &in, //When you use an array as an argument to the top-level function. V. HLS assumes Memory is off-chip
		hls::stream <mat> &out//Streams must be passed in and out of functions by-reference
			)
{
#pragma HLS ARRAY_PARTITIOn dim=1 type=complete variable=coeffs_2D
	int index;
	int counter;
	int cnt=0;
	mat input;
	mat output;
/*#ifndef __SYNTHESIS__
	FIFO_Init:for(int i=0; i<STREAM_DEPTH;i++){
		fifo_0.write((mat)0);
		fifo_1.write((mat)0);
		fifo_2.write((mat)0);
		fifo_3.write((mat)0);
	}
//#endif*/
	Height:for (int i=0;i<HEIGHT+2;i++){
		Width:for (int j=0;j<WIDTH;j+=I){
#pragma HLS PIPELINE II=1
			//Spatial_Parallelism:for (int v=0;v<I;v++){ //It is not actually a Parallelism since one result must come after the previous
			if (i<HEIGHT){
				in.read(input);
			}
			else{input=0;}
			if(i==0 && j==0){
				index=0;
			}
			else if(i==0 && j==1){
				index=1;
			}
			else if(i==0 && j==2){
				index=2;
			}
			else if(i==0 && j>2 && j<WIDTH-3){
				index=3;
			}
			else if(i==0 && j==WIDTH-3){
				index=4;
			}
			else if(i==0 && j==WIDTH-2){
				index=5;
			}
			else if(i==0 && j==WIDTH-1){
				index=6;
			}
			else if(i==1 && j==0){
				index=7;
			}
			else if(i==1 && j==1){
				index=8;
			}
			else if(i==1 && j==2){
				index=9;
			}
			else if(i==1 && j>2 && j<WIDTH-3){
				index=10;
			}
			else if(i==1 && j==WIDTH-3){
				index=11;
			}
			else if(i==1 && j==WIDTH-2){
				index=12;
			}
			else if(i==1 && j==WIDTH-1){
				index=13;
			}
			else if(i==2 && j==0){
				index=14;
			}
			else if(i==2 && j==1){
				index=15;
			}
			else if(i==2 && j==2){
				index=16;
			}
			else if(i==2 && j>2 && j<WIDTH-3){
				index=17;
			}
			else if(i==2 && j==WIDTH-3){
				index=18;
			}
			else if(i==2 && j==WIDTH-2){
				index=19;
			}
			else if(i==2 && j==WIDTH-1){
				index=20;
			}
			else if(i>2 && i<HEIGHT-3 && j==0){
				index=21;
			}
			else if(i>2 && i<HEIGHT-3 && j==1){
				index=22;
			}
			else if(i>2 && i<HEIGHT-3 && j==2){
				index=23;
			}
			else if(i>2 && i<HEIGHT-3 && j==WIDTH-3){
				index=25;
			}
			else if(i>2 && i<HEIGHT-3 && j==WIDTH-2){
				index=26;
			}
			else if(i>2 && i<HEIGHT-3 && j==WIDTH-1){
				index=27;
			}
			else if(i==HEIGHT-1 && j==0){
				index=42;
			}
			else if(i==HEIGHT-1 && j==1){
				index=43;
			}
			else if(i==HEIGHT-1 && j==2){
				index=44;
			}
			else if(i==HEIGHT-2 && j==0){
				index=35;
			}
			else if(i==HEIGHT-2 && j==1){
				index=36;
			}
			else if(i==HEIGHT-2 && j==2){
				index=37;
			}
			else if(i==HEIGHT-3 && j==0){
				index=28;
			}
			else if(i==HEIGHT-3 && j==1){
				index=29;
			}
			else if(i==HEIGHT-3 && j==2){
				index=30;
			}
			else if(i==HEIGHT-1 && j>2 && j<WIDTH-3){
				index=45;
			}
			else if(i==HEIGHT-2 && j>2 && j<WIDTH-3){
				index=38;
			}
			else if(i==HEIGHT-3 && j>2 && j<WIDTH-3){
				index=31;
			}
			else if(i==HEIGHT-1 && j==WIDTH-1){
				index=48;
			}
			else if(i==HEIGHT-1 && j==WIDTH-2){
				index=47;
			}
			else if(i==HEIGHT-1 && j==WIDTH-3){
				index=46;
			}
			else if(i==HEIGHT-2 && j==WIDTH-1){
				index=41;
			}
			else if(i==HEIGHT-2 && j==WIDTH-2){
				index=40;
			}
			else if(i==HEIGHT-2 && j==WIDTH-3){
				index=39;
			}
			else if(i==HEIGHT-3 && j==WIDTH-1){
				index=34;
			}
			else if(i==HEIGHT-3 && j==WIDTH-2){
				index=33;
			}
			else if(i==HEIGHT-3 && j==WIDTH-3){
				index=32;
			}
			else {  //if((i>2 && i<HEIGHT-3 && j>2 && j<WIDTH-3)||(i>HEIGHT-1)){
				index=24;
			}
			DCMI_Compute_Engine_Core(input,output,index);
			if (i*WIDTH+j>2*WIDTH+1){ //D=2 so i-D and j+v-D would also be correct
				out.write(output);
				/*if(j<2){
					out[i-3][WIDTH-2+(j)]=output;
				}
				else{
				out[i-2][j-2]=output;
				}*/

			}
		}
	}
/*
	DCMI_Compute_Engine_Core(input,output,index,(HEIGHT+2)*WIDTH);
	out.write(output);
	DCMI_Compute_Engine_Core(input,output,index,(HEIGHT+2)*WIDTH+1);
	out.write(output);
*/
	/*//input[I]={0.0};
	DCMI_Compute_Engine_Core(HEIGHT+1,WIDTH,input,output,coeffs);
	out[HEIGHT-1][WIDTH-2]=output[0];
	DCMI_Compute_Engine_Core(HEIGHT+1,WIDTH+1,input,output,coeffs);
	out[HEIGHT-1][WIDTH-1]=output[0];*/
}
