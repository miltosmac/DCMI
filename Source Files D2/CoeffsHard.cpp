#include "DCMI_Jacobi.h"

mat coef0_0[N] = { 0.0f,0.0f,0.0f,0.0f,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0.123456791f,0.0123456791f,0,0,0,0.0123456791f,0.0123456791f};
mat coef0_1[N] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0.135802478f,0.135802478f,0.0123456791f,0,0,0.0246913582f,0.0246913582f,0.0123456791f };
mat coef0_2[N] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0.135802478f,0.148148149f,0.135802478f,0.0123456791f,0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f};
mat coef0_3W4[N] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0.0123456791f,0.135802478f,0.148148149f,0.135802478f,0.0123456791f,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f };
mat coef1_0[N] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,0.135802478f,0.0246913582f,0,0,0,0.135802478f,0.0246913582f,0,0,0,0.0123456791f,0.0123456791f };
mat coef1_1[N] = { 0,0,0,0,0,0,0,0,0,0,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0246913582f,0.0246913582f,0.0123456791f };
mat coef1_2[N] = { 0,0,0,0,0,0,0,0,0,0,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f };
mat coef1_3W4[N] = { 0,0,0,0,0,0,0,0,0,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f };
mat coef2_0[N] = { 0,0,0,0,0,0,0,0,0.135802478f,0.0246913582f,0,0,1,0.148148149f,0.0370370373f,0,0,0,0.135802478f,0.0246913582f,0,0,0,0.0123456791f,0.0123456791f };
mat coef2_1[N] = { 0,0,0,0,0,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0740740746f,0.0740740746f,0.0370370373f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0246913582f,0.0246913582f,0.0123456791f };
mat coef2_2[N] = { 0,0,0,0,0,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0740740746f,0.111111112f,0.0740740746f,0.0370370373f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f};
mat coef2_3W4[N] = { 0,0,0,0,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0370370373f,0.0740740746f,0.111111112f,0.0740740746f,0.0370370373f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f };
mat coef3H4_0[N] = { 0,0,0,0.0123456791f,0.0123456791f,0,0,0,0.135802478f,0.0246913582f,0,0,1,0.148148149f,0.0370370373f,0,0,0,0.135802478f,0.0246913582f,0,0,0,0.0123456791f,0.0123456791f };
mat coef3H4_1[N] = { 0,0,0.0246913582f,0.0246913582f,0.0123456791f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0740740746f,0.0740740746f,0.0370370373f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0246913582f,0.0246913582f,0.0123456791f};
mat coef3H4_2[N] = { 0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0740740746f,0.111111112f,0.0740740746f,0.0370370373f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f };
mat coeffs[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0370370373f,0.0740740746f,0.111111112f,0.0740740746f,0.0370370373f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f };

mat coef0_W1 [N]= { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0.0123456791f,0.123456791f,0,0,0,0.0123456791f,0.0123456791f,0,0,0};
mat coef0_W2 [N]= { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0.0123456791f,0.135802478f,0.135802478f,0,0,0.0123456791f,0.0246913582f,0.0246913582f,0,0 };
mat coef0_W3[N] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0.0123456791f,0.135802478f,0.148148149f,0.135802478f,0,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0};
mat coef1_W1[N] = { 0,0,0,0,0,0,0,0,0,0,0.0246913582f,0.135802478f,1,0,0,0.0246913582f,0.135802478f,0,0,0,0.0123456791f,0.0123456791f,0,0,0 };
mat coef1_W2[N] = { 0,0,0,0,0,0,0,0,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0123456791f,0.0246913582f,0.0246913582f,0,0 };
mat coef1_W3[N] = { 0,0,0,0,0,0,0,0,0,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0 };
mat coef2_W1[N] = { 0,0,0,0,0,0.0246913582f,0.135802478f,0,0,0,0.0370370373f,0.148148149f,1,0,0,0.0246913582f,0.135802478f,0,0,0,0.0123456791f,0.0123456791f,0,0,0 };
mat coef2_W2 [N]= { 0,0,0,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0370370373f,0.0740740746f,0.0740740746f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0123456791f,0.0246913582f,0.0246913582f,0,0 };
mat coef2_W3 [N]= { 0,0,0,0,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0370370373f,0.0740740746f,0.111111112f,0.0740740746f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0 };
mat coef3H4_W1 [N]= { 0.0123456791f,0.0123456791f,0,0,0,0.0246913582f,0.135802478f,0,0,0,0.0370370373f,0.148148149f,1,0,0,0.0246913582f,0.135802478f,0,0,0,0.0123456791f,0.0123456791f,0,0,0};
mat coef3H4_W2 [N]= { 0.0123456791f,0.0246913582f,0.0246913582f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0370370373f,0.0740740746f,0.0740740746f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0123456791f,0.0246913582f,0.0246913582f,0,0 };
mat coef3H4_W3 [N]= { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0370370373f,0.0740740746f,0.111111112f,0.0740740746f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0 };

mat coefH1_0[N] = { 0,0,0,0.0123456791f,0.0123456791f,0,0,0,0.123456791f,0.0123456791f,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
mat coefH1_1[N] = { 0,0,0.0246913582f,0.0246913582f,0.0123456791f,0,0,0.135802478f,0.135802478f,0.0123456791f,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
mat coefH1_2[N] = { 0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0,0.135802478f,0.148148149f,0.135802478f,0.0123456791f,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH1_3W4[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0.0123456791,0.135802478f,0.148148149f,0.135802478f,0.0123456791f,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_0[N] = { 0,0,0,0.0123456791f,0.0123456791f,0,0,0,0.135802478f,0.0246913582f,0,0,1,0.135802478f,0.0246913582f,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_1[N] = { 0,0,0.0246913582f,0.0246913582f,0.0123456791f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_2[N] = { 0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_3W4[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0,0,0,0,0,0,0,0,0 };
mat coefH3_0[N] = { 0,0,0,0.0123456791f,0.0123456791f,0,0,0,0.135802478f,0.0246913582f,0,0,1,0.148148149f,0.0370370373f,0,0,0,0.135802478f,0.0246913582f,0,0,0,0,0 };
mat coefH3_1[N] = { 0,0,0.0246913582f,0.0246913582f,0.0123456791f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0.0740740746f,0.0740740746f,0.0370370373f,0,0,0.0493827164f,0.0493827164f,0.0246913582f,0,0,0,0,0 };
mat coefH3_2[N] = { 0,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0.0740740746f,0.111111112f,0.0740740746f,0.0370370373f,0,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0,0,0,0 };
mat coefH3_3W4[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0.0123456791f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0.0370370373f,0.0740740746f,0.111111112f,0.0740740746f,0.0370370373f,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0.0246913582f,0,0,0,0,0 };

mat coefH1_W1[N] = { 0.0123456791f,0.0123456791f,0,0,0,0.0123456791f,0.123456791f,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH1_W2[N] = { 0.0123456791f,0.0246913582f,0.0246913582f,0,0,0.0123456791f,0.135802478f,0.135802478f,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH1_W3[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0,0.0123456791f,0.135802478f,0.148148149f,0.135802478f,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_W1[N] = { 0.0123456791f,0.0123456791f,0,0,0,0.0246913582f,0.135802478f,0,0,0,0.0246913582f,0.135802478f,1,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_W2[N] = { 0.0123456791f,0.0246913582f,0.0246913582f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH2_W3[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0,0,0,0,0,0,0,0,0,0 };
mat coefH3_W1[N] = { 0.0123456791f,0.0123456791f,0,0,0,0.0246913582f,0.135802478f,0,0,0,0.0370370373f,0.148148149f,1,0,0,0.0246913582f,0.135802478f,0,0,0,0,0,0,0,0 };
mat coefH3_W2[N] = { 0.0123456791f,0.0246913582f,0.0246913582f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0.0370370373f,0.0740740746f,0.0740740746f,0,0,0.0246913582f,0.0493827164f,0.0493827164f,0,0,0,0,0,0,0 };
mat coefH3_W3[N] = { 0.0123456791f,0.0246913582f,0.0370370373f,0.0246913582f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0.0370370373f,0.0740740746f,0.111111112f,0.0740740746f,0,0.0246913582f,0.0493827164f,0.0740740746f,0.0493827164f,0,0,0,0,0,0 };
