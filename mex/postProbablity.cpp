#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"


#define IN_x        prhs[0]
#define IN_y        prhs[1]
#define IN_sigma2   prhs[2]
#define IN_outlier  prhs[3]
#define OUT_P       plhs[0]
#define PI          3.14159265358979323846

void displ(double *P ,int N, int M){
    for(int m = 0; m < M; ++m){
        for(int n = 0; n < N; ++n){
            printf("%.8f", *(P + m + n * M));
            if(n < N - 1) {
                printf(" ");
            }
        }
        printf("\n");
    }
}

void cpd_comp(double *x, double *y, double *sigma2, double *outliers, double * P, int N ,int M, int D){
    double k_sig = -2 * (*sigma2);
    double outlier_tmp = ((*outliers) * M * pow(-k_sig * PI, 0.5 * D)) / (N * (1 - (*outliers)));

    for(int n = 0; n < N; ++n){
        double sp = 0;
        for(int m = 0; m < M; ++m){
            /*printf("----------------------------------\n");
            printf("n:%d m:%d\n",n + 1,m + 1);*/
            double diff = 0;
            double razn = 0;
            for(int d = 0; d < D; ++d){
                diff = *(x + n + d * N) - *(y + m + d * M);
                diff = diff * diff;
                razn += diff;
              /*  printf("diff : %.2f\n", diff);
                printf("len : %.2f\n", len);*/

            }
            *(P + m + n * M) = exp(razn / k_sig);
            sp += *(P + m + n * M);
           /* printf("P_(%d)(%d) : %.8f\n", m + 1, n + 1, *(P + m + n * M));
            printf("res : %.8f\n", res);*/
        }
        sp += outlier_tmp;
       // printf("%d : %f\n", n + 1, sp);
        for(int m = 0; m < M; ++m){
            *(P + m + n * M) = *(P + m + n * M) / sp;
        }
//        printf("----------------------------------\n");
    }
    //displ(P, N, M);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *x, *y, *sigma2, *outliers, *P;
    int N, M, D;
    // input
    N = (int) mxGetM(IN_x);
    M = (int) mxGetM(IN_y);
    D = (int) mxGetN(IN_x);
    x = mxGetPr(IN_x);
    y = mxGetPr(IN_y);
    sigma2 = mxGetPr(IN_sigma2);
    outliers = mxGetPr(IN_outlier);

    //output
    OUT_P = mxCreateDoubleMatrix((mwSize) M, (mwSize) N, mxREAL);
    P = mxGetPr(OUT_P);
    cpd_comp(x, y, sigma2, outliers, P, N, M, D);
}
