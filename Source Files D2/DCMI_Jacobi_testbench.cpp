/***********************************
9-Point Jacobi 2-D  DCMI TestBench File
Miliadis Mataragkas, February 2022
University of Patras
***********************************/
#include <iostream>
#include <math.h>
#include <fstream>
#include "DCMI_Jacobi.h"

using namespace std;

int main(int argc, char **argv)
{
	mat temp;
	mat c=1.0f/9.0f;
	int err_cnt = 0;
	float mult;
	mat A[HEIGHT][WIDTH],sw_result[HEIGHT][WIDTH],hw_result[HEIGHT][WIDTH];
	mat coeff[2*R*D+1];
	hls::stream<mat> testream_in ("instream");
	hls::stream<mat> testream_out ("outstream");
	hls::stream<mat> coefficient_in ("load_coeffs");
	ofstream txtfile;
	txtfile.open ("input_mat.txt");

	for (int i=0; i<HEIGHT;i++){
		for(int j=0; j<WIDTH; j++){
			mult=float(rand());
			//A[i][j]= mult*float(rand())/float(RAND_MAX); //Initialize A matrix as random floats
			A[i][j]=static_cast <mat>(i*WIDTH+j);
			testream_in.write(A[i][j]);
			txtfile << A[i][j]<<"  ";
		}
		txtfile <<"\n";
	}
	txtfile.close();
	//Calculate The Appropriate Coefficient Arrays for D=2

	//Coefficient_Generation();

	// Calculate DCMI Result

	jacobi9d(testream_in,testream_out);

	//Generate the expected results
	//Iterate over the rows of the 2D input matrix
	for(int t=0;t<D;t++){
		for (int i=0; i<HEIGHT; i++){
			//Iterate the columns of the 2D input matrix
			for (int j=0; j<WIDTH; j++){
				if (i>0 && i < HEIGHT-1 && j>0 && j < WIDTH-1){
					temp= A[i-1][j-1]+A[i-1][j]+A[i-1][j+1]+ // First Line
							A[i][j-1]+A[i][j]+A[i][j+1]+      // Second Line
							A[i+1][j-1]+A[i+1][j]+A[i+1][j+1];   // Third Line
					sw_result[i][j]=c*temp; //Create SW result
				}else {
					sw_result[i][j]=A[i][j]; //Results for the HALO
				}
			}
		}
		for (int i=0;i<HEIGHT;i++){
			for (int j=0;j<WIDTH;j++){
				A[i][j]=sw_result[i][j];
			}
		}
	}
	/* The coefficient calculation is carried out at design time
	 *  using a regular computer so its performance is not a significant concern.*/
	txtfile.open ("output_mat.txt");
	for (int i=0; i<HEIGHT; i++){
		for (int j=0; j<WIDTH; j++){
		// Check HW result against SW
				if (!testream_out.empty()){testream_out.read(hw_result[i][j]);} //Values read from stream and into the HW result variable for checking
				txtfile << hw_result[i][j]<<"  ";
			if (fabs(float(hw_result[i][j] - sw_result[i][j]))>1e-0f) { //Just checking the Kernel results not the full size matrix with the HALO
 				cout << "At iteration i:" << i << " and j:" <<j << endl;
				err_cnt++;
			}
		}
		txtfile <<"\n";
	}
	txtfile.close();
	if (err_cnt){
		cout << "ERROR: " << err_cnt << " mismatches detected!" << endl;
	}else{
		cout << "Test passes." << endl;
	}
}
