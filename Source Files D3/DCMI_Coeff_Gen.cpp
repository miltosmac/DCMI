#include "DCMI_Jacobi.h"
#include <fstream>
#include <iostream>
#include <string>


void Coefficient_Generation(
	)
{
	int div_i, mod_i, counter;
	const int L = 2 * R * D + 1; //Length of array of affected values
	const mat c = 1.0f / 9.0f;
	int div, rem,temp_i,temp_j;
	mat  T[2 * R * D + 1] = { 0 };
	auto top_left = new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto top_right= new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto bot_left = new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto bot_right= new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto c_temp   = new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto bl_temp  = new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto tr_temp  = new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	auto br_temp  = new mat[(D + 2) * (D + 2)][2 * R * D + 1];
	std::ofstream txtfile;
	std::ofstream coeffs_txt;
	
	//Initialize the Coeff Array
	for (int i = 0; i < ((D + 2) * (D + 2)); i++) {
		for (int j = 0; j < L; j++) {
			if (j == R * D) {
				bot_left[i][j] = 1.0f;
				bot_right[i][j] = 1.0f;
				top_right[i][j] = 1.0f;
				top_left[i][j] = 1.0f;
			}
			else {
				bot_left[i][j] = 0.0f;
				bot_right[i][j] = 0.0f;
				top_left[i][j] = 0.0f;
				top_right[i][j] = 0.0f;
			}
		}
	}
	coeffs_txt.open("coeffs_file.txt");
	/*Generate the Coefficient Array\
	* for the TOP-LEFT corner */
	for (int t = 0; t < D; t++) {
		for (int ii = 0; ii <= t + 2; ii++) {
			for (int jj = 0; jj <= t + 2; jj++) {
				if ((ii == 0 || jj == 0) && (t == D - 1)) {
					T[R * D] = 1.0f;
				}
				for (int i = 0; i < L; i++) {
					for (int j = -1; j < 2; j++) {
						/*Keep in mind that that the data element accessed corresponds to RB[RD], the top left element of the stencil pattern corresponds to RB[0] and the last to RB[2RD]*/
						div = i / WIDTH; //i
						rem = i % WIDTH; //j
						if (ii == t + 2)
						{
							temp_i = ii - 1;
						}
						else
						{
							temp_i = ii;
						}
						if (jj == t + 2)
						{
							temp_j = jj - 1;
						}
						else
						{
							temp_j = jj;
						}
						if (div > D - ii && rem > D - jj) {
							if ((i + j - WIDTH) > -1 && (i + j - WIDTH) < L) {
								T[i] += c * top_left[temp_i * (2 + D) + temp_j][i + j - WIDTH];
							}
							if ((i + j) > -1 && (i + j) < L) {
								T[i] += c * top_left[temp_i * (2 + D) + temp_j][i + j];
							}
							if ((i + j + WIDTH) > -1 && (i + j + WIDTH) < L) {
								T[i] += c * top_left[temp_i * (2 + D) + temp_j][i + j + WIDTH];
							}
							if (t == D - 1 && ((ii == 0 || jj == 0) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
							if (t == D - 2 && ((ii == 0 || jj == 0) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c ;
							}
						}
					}
				}
				for (int i = 0; i < L; i++) {
					c_temp[ii * (2 + D) + jj][i] = T[i];
					T[i] = 0;
				}
			}
		}
		for (int ii = 0; ii <= t + 2; ii++) {
			for (int jj = 0; jj <= t + 2; jj++) {
				txtfile.open("i" + std::to_string(ii) + "j" + std::to_string(jj) + ".txt");
				if (t == D - 1) { coeffs_txt << "{"; }
				for (int i = 0; i < L; i++) {
					top_left[ii * (2 + D) + jj][i] = c_temp[ii * (2 + D) + jj][i];
					mod_i = i % WIDTH;
					if (t == D - 1 &&mod_i < (D + 4)) {
						txtfile << i << ":" << top_left[ii * (2 + D) + jj][i] << "\n";
						coeffs_txt << std::setprecision(9) << top_left[ii * (2 + D) + jj][i] << ",";
					}
				}
				if (t == D - 1) { coeffs_txt << "}\n"; }
				txtfile.close();
			}
		}
	}
	for (int i = 0; i < L; i++) {
		div_i = i / WIDTH;
		mod_i = i % WIDTH;
		counter = div_i * (D + 4) + mod_i;
		if (counter < N) {
			coef0_0[counter] = top_left[0][i];
			coef0_1[counter] = top_left[1][i];
			coef0_2[counter] = top_left[2][i];
			coef0_3[counter] = top_left[3][i];
			coef0_4W5[counter] = top_left[4][i];
			coef1_0[counter] = top_left[5][i];
			coef1_1[counter] = top_left[6][i];
			coef1_2[counter] = top_left[7][i];
			coef1_3[counter] = top_left[8][i];
			coef1_4W5[counter] = top_left[9][i];
			coef2_0[counter] = top_left[10][i];
			coef2_1[counter] = top_left[11][i];
			coef2_2[counter] = top_left[12][i];
			coef2_3[counter] = top_left[13][i];
			coef2_4W5[counter] = top_left[14][i];
			coef3_0[counter] = top_left[15][i];
			coef3_1[counter] = top_left[16][i];
			coef3_2[counter] = top_left[17][i];
			coef3_3[counter] = top_left[18][i];
			coef3_4W5[counter] = top_left[19][i];
			coef4H5_0[counter] = top_left[20][i];
			coef4H5_1[counter] = top_left[21][i];
			coef4H5_2[counter] = top_left[22][i];
			coef4H5_3[counter] = top_left[23][i];
			coeffs[counter] = top_left[24][i];
		}
	}
	delete[] c_temp;
	delete[] top_left;
	/*Generate the Coefficient Array\
		* for the TOP-RIGHT corner */
	for (int t = 0; t < D; t++) {
		for (int ii = 0; ii <= t+2; ii++) {
			for (int jj = WIDTH-1; jj > WIDTH-4-t; jj--) {
				if ((ii == 0 || jj == WIDTH-1) && (t == D-1)){
					T[R * D] = 1.0f;
				}
				for (int i = 0; i < L; i++) {
					for (int j = -1; j < 2; j++) {
						/*Keep in mind that that the data element accessed corresponds to RB[R*D], the top left element of the stencil pattern corresponds to RB[0] and the last to RB[2RD]*/
						div = i / WIDTH; //i
						rem = i % WIDTH; //j
						if (ii == t + 2){
							temp_i = ii - 1;
						}
						else{
							temp_i = ii;
						}
						if (jj == WIDTH - 3 - t){
							temp_j = WIDTH - 1 - jj - 1;
						}
						else{
							temp_j = WIDTH - 1- jj;
						}
						if (div > D - ii && rem < (WIDTH - 1) + D - jj) {
							if ((i + j - WIDTH) > -1 && (i + j - WIDTH) < L) {
								T[i] += c * top_right[temp_i * (2 + D) + temp_j][i + j - WIDTH];
							}
							if ((i + j) > -1 && (i + j) < L) {
								T[i] += c * top_right[temp_i * (2 + D) + temp_j][i + j];
							}
							if ((i + j + WIDTH) > -1 && (i + j + WIDTH) < L) {
								T[i] += c * top_right[temp_i * (2 + D) + temp_j][i + j + WIDTH];
							}
							if (t == D - 1 && ((ii == 0 || jj == WIDTH-1) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
							if (t == D - 2 && ((ii == 0 || jj == WIDTH - 1) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
						}
					}
				}
				for (int i = 0; i < L; i++) {
					tr_temp[ii*(2+D)+ WIDTH - 1 - jj][i] = T[i];
					T[i] = 0;
				}
			}
		}
		for (int ii = 0; ii <= t + 2; ii++) {
			for (int jj = WIDTH - 1; jj > WIDTH - 4 - t; jj--) {
				txtfile.open("i" + std::to_string(ii) + "j" + std::to_string(jj) + ".txt");
				if (t == D - 1) { coeffs_txt << "{"; }
				for (int i = 0; i < L; i++) {
					top_right[ii * (2 + D) + WIDTH-1-jj][i] = tr_temp[ii * (2 + D) + WIDTH - 1 - jj][i];
					mod_i = i % WIDTH;
					if (t == D - 1 && mod_i <(D+4)) {
						txtfile << i << ":" << top_right[ii * (2 + D) + WIDTH - 1 - jj][i] << "\n";
						coeffs_txt << std::setprecision(9) << top_right[ii * (2 + D) + WIDTH - 1 - jj][i] << ",";
					}
				}
				txtfile.close();
				if (t == D - 1) { coeffs_txt << "}\n"; }
			}
		}
		coeffs_txt << "\n";
	}
	for (int i = 0; i < L; i++) {
		div_i = i / WIDTH;
		mod_i = i % WIDTH;
		counter = div_i * (D + 4) + mod_i;
		if (counter < N) {
			coef0_W1[counter] = top_right[0][i];
			coef0_W2[counter] = top_right[1][i];
			coef0_W3[counter] = top_right[2][i];
			coef0_W4[counter] = top_right[3][i];
			coef1_W1[counter] = top_right[5][i];
			coef1_W2[counter] = top_right[6][i];
			coef1_W3[counter] = top_right[7][i];
			coef1_W4[counter] = top_right[8][i];
			coef2_W1[counter] = top_right[10][i];
			coef2_W2[counter] = top_right[11][i];
			coef2_W3[counter] = top_right[12][i];
			coef2_W4[counter] = top_right[13][i];
			coef3_W1[counter] = top_right[15][i];
			coef3_W2[counter] = top_right[16][i];
			coef3_W3[counter] = top_right[17][i];
			coef3_W4[counter] = top_right[18][i];
			coef4H5_W1[counter] = top_right[20][i];
			coef4H5_W2[counter] = top_right[21][i];
			coef4H5_W3[counter] = top_right[22][i];
			coef4H5_W4[counter] = top_right[23][i];
		}
	}
	delete[] top_right;
	delete[] tr_temp;

	/*Generate the Coefficient Array\
		* for the BOTTOM-LEFT corner */
	for (int t = 0; t < D; t++) {
		for (int ii = HEIGHT -1 ; ii  >= HEIGHT-3-t; ii--) {
			for (int jj = 0; jj <=t+2 ; jj++) {
				if ((ii == HEIGHT-1 || jj == 0) && (t== D-1)) {
					T[R * D] = 1.0f;
				}
				for (int i = 0; i < L; i++) {
					for (int j = -1; j < 2; j++) {
						/*Keep in mind that that the data element accessed corresponds to RB[RD], the top left element of the stencil pattern corresponds to RB[0] and the last to RB[2RD]*/
						div = i / WIDTH; //i
						rem = i % WIDTH; //j
						if (ii == (HEIGHT-3-t)){
							temp_i = HEIGHT-1- ii -1;
						}
						else{
							temp_i = HEIGHT-1-ii;
						}
						if (jj == t + 2){
							temp_j = jj - 1;
						}
						else{
							temp_j = jj;
						}
						if (div < (HEIGHT-1) + D - ii && rem > D - jj) {
							if ((i + j - WIDTH) > -1 && (i + j - WIDTH) < L) {
								T[i] += c * bot_left[temp_i * (2 + D) + temp_j][i + j - WIDTH];
							}
							if ((i + j) > -1 && (i + j) < L) {
								T[i] += c * bot_left[temp_i * (2 + D) + temp_j][i + j];
							}
							if ((i + j + WIDTH) > -1 && (i + j + WIDTH) < L) {
								T[i] += c * bot_left[temp_i * (2 + D) + temp_j][i + j + WIDTH];
							}
							if (t == D - 1 && ((ii == HEIGHT-1 || jj == 0) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
							if (t == D - 2 && ((ii == HEIGHT - 1  || jj == 0) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
						}
					}
				}
				for (int i = 0; i < L; i++) {
					bl_temp[(HEIGHT-1-ii)*(2+D)+jj][i] = T[i];
					T[i] = 0;
				}
			}
		}
		for (int ii = HEIGHT - 1; ii > HEIGHT - 4 - t; ii--) {
			for (int jj = 0; jj <= t + 2; jj++) {
				if (t == D - 1) { coeffs_txt << "{"; }
				txtfile.open("i" + std::to_string(ii) + "j" + std::to_string(jj) + ".txt");
				for (int i = 0; i < L; i++) {
					bot_left[(HEIGHT - 1 - ii) * (2 + D) + jj][i] = bl_temp[(HEIGHT - 1 - ii) * (2 + D) + jj][i];
					mod_i = i % WIDTH;
					if (t == D - 1 && mod_i <(D+4)) {
						txtfile << i << ":" << bot_left[(HEIGHT - 1 - ii) * (2 + D) + jj][i] << "\n";
						coeffs_txt << std::setprecision(9) << bot_left[(HEIGHT - 1 - ii) * (2 + D) + jj][i] << ",";
					}
				}
				txtfile.close();
				if (t == D - 1) { coeffs_txt << "}\n"; }
			}
		}
		coeffs_txt << "\n";
	}
	for (int i = 0; i < L; i++) {
		div_i = i / WIDTH;
		mod_i = i % WIDTH;
		counter = div_i * (D + 4) + mod_i;
		if (counter < N) {
			coefH1_0[counter] = bot_left[0][i];
			coefH1_1[counter] = bot_left[1][i];
			coefH1_2[counter] = bot_left[2][i];
			coefH1_3[counter] = bot_left[3][i];
			coefH1_4W5[counter] = bot_left[4][i];
			coefH2_0[counter] = bot_left[5][i];
			coefH2_1[counter] = bot_left[6][i];
			coefH2_2[counter] = bot_left[7][i];
			coefH2_3[counter] = bot_left[8][i];
			coefH2_4W5[counter] = bot_left[9][i];
			coefH3_0[counter] = bot_left[10][i];
			coefH3_1[counter] = bot_left[11][i];
			coefH3_2[counter] = bot_left[12][i];
			coefH3_3[counter] = bot_left[13][i];
			coefH3_4W5[counter] = bot_left[14][i];
			coefH4_0[counter] = bot_left[15][i];
			coefH4_1[counter] = bot_left[16][i];
			coefH4_2[counter] = bot_left[17][i];
			coefH4_3[counter] = bot_left[18][i];
			coefH4_4W5[counter] = bot_left[19][i];
		}
	}
	delete[] bot_left;
	delete[] bl_temp;
	/*Generate the Coefficient Array\
	* for the BOTTOM-RIGHT corner */
	for (int t = 0; t < D; t++) {
		for (int ii = HEIGHT - 1; ii >= HEIGHT - 3 - t; ii--) {
			for (int jj = WIDTH - 1; jj > WIDTH - 4 - t; jj--) {
				if ((ii == HEIGHT - 1 || jj == WIDTH - 1) && (t == D - 1)) {
					T[R * D] = 1.0f;
				}
				for (int i = 0; i < L; i++) {
					for (int j = -1; j < 2; j++) {
						/*Keep in mind that that the data element accessed corresponds to RB[RD], the top left element of the stencil pattern corresponds to RB[0] and the last to RB[2RD]*/
						div = i / WIDTH; //i
						rem = i % WIDTH; //j
						if (ii == (HEIGHT - 3 - t)) {
							temp_i = HEIGHT - 1 - ii - 1;
						}
						else {
							temp_i = HEIGHT - 1 - ii;
						}
						if (jj == (WIDTH - 3 - t)) {
							temp_j = WIDTH -1 - jj - 1;
						}
						else {
							temp_j = WIDTH - 1 - jj;
						}
						if (div < (HEIGHT - 1) + D - ii && rem < (WIDTH - 1) + D - jj) {
							if ((i + j - WIDTH) > -1 && (i + j - WIDTH) < L) {
								T[i] += c * bot_right[temp_i * (2 + D) + temp_j][i + j - WIDTH];
							}
							if ((i + j) > -1 && (i + j) < L) {
								T[i] += c * bot_right[temp_i * (2 + D) + temp_j][i + j];
							}
							if ((i + j + WIDTH) > -1 && (i + j + WIDTH) < L) {
								T[i] += c * bot_right[temp_i * (2 + D) + temp_j][i + j + WIDTH];
							}
							if (t == D - 1 && ((ii == HEIGHT - 1 || jj == WIDTH - 1) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
							if (t == D - 2 && ((ii == HEIGHT - 1 || jj == WIDTH - 1) && abs(div - D) <= 1 && abs(rem - D) <= 1) && j == 0) {
								T[i] += c;
							}
						}
					}
				}
				for (int i = 0; i < L; i++) {
					br_temp[(HEIGHT - 1 - ii) * (2 + D) + WIDTH-1-jj][i] = T[i];
					T[i] = 0;
				}
			}
		}
		for (int ii = HEIGHT - 1; ii > HEIGHT - 4 - t; ii--) {
			for (int jj = WIDTH - 1; jj > WIDTH - 4 - t; jj--) {
				if (t == D - 1) { coeffs_txt << "{"; }
				txtfile.open("i" + std::to_string(ii) + "j" + std::to_string(jj) + ".txt");
				for (int i = 0; i < L; i++) {
					mod_i = i % WIDTH;
					bot_right[(HEIGHT - 1 - ii) * (2 + D) + WIDTH - 1 - jj][i] = br_temp[(HEIGHT - 1 - ii) * (2 + D) + WIDTH - 1 - jj][i];
					if (t == D - 1 && mod_i <(D+4)) {
						txtfile << i << ":" << bot_right[(HEIGHT - 1 - ii) * (2 + D) + WIDTH - 1 - jj][i] << "\n";
						coeffs_txt << std::setprecision(9) << bot_right[(HEIGHT - 1 - ii) * (2 + D) + WIDTH - 1 - jj][i] << ",";
					}
				}
				txtfile.close();
				if (t == D - 1) { coeffs_txt << "}\n"; }
			}
		}
	}
	for (int i = 0; i < L; i++) {
		div_i = i / WIDTH;
		mod_i = i % WIDTH;
		counter = div_i * (D + 4) + mod_i;
		if (counter < N) {
			coefH1_W1[counter] = bot_right[0][i];
			coefH1_W2[counter] = bot_right[1][i];
			coefH1_W3[counter] = bot_right[2][i];
			coefH1_W4[counter] = bot_right[3][i];
			coefH2_W1[counter] = bot_right[5][i];
			coefH2_W2[counter] = bot_right[6][i];
			coefH2_W3[counter] = bot_right[7][i];
			coefH2_W4[counter] = bot_right[8][i];
			coefH3_W1[counter] = bot_right[10][i];
			coefH3_W2[counter] = bot_right[11][i];
			coefH3_W3[counter] = bot_right[12][i];
			coefH3_W4[counter] = bot_right[13][i];
			coefH4_W1[counter] = bot_right[15][i];
			coefH4_W2[counter] = bot_right[16][i];
			coefH4_W3[counter] = bot_right[17][i];
			coefH4_W4[counter] = bot_right[18][i];
		}
	}
	delete[] bot_right;
	delete[] br_temp;
	coeffs_txt.close();
}
