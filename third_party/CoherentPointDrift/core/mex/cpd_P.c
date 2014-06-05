/*
Andriy Myronenko
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))

void cpd_comp(
		double* x,
		double* y, 
        double* sigma2,
		double* outlier,
        double* P1,
        double* Pt1,
        double* Px,
    	double* E,
        int N,
		int M,
        int D
        )

{
  int		n, m, d;
  double	ksig, diff, razn, outlier_tmp, sp;
  double	*P, *temp_x;
  
  P = (double*) calloc(M, sizeof(double));
  temp_x = (double*) calloc(D, sizeof(double));
  
  ksig = -2.0 * *sigma2;
  outlier_tmp=(*outlier*M*pow (-ksig*3.14159265358979,0.5*D))/((1-*outlier)*N); // right half of the denominator
 /* printf ("ksig = %lf\n", *sigma2);*/
  /* outlier_tmp=*outlier*N/(1- *outlier)/M*(-ksig*3.14159265358979); */
  
  // for each point in the target set
  for (n=0; n < N; n++) {
      
      sp=0;
      // for each point in the input set
      for (m=0; m < M; m++) {
          razn=0;
          
          // for each dimension (x, y, z)
          for (d=0; d < D; d++) {
             // diff is |Xn - Ym|^2
             diff=*(x+n+d*N)-*(y+m+d*M);  diff=diff*diff;
             // razn is 
             razn+=diff;
          }
          
          // sp = P[m] = e^(|Xn - Ym|^2 / (-2.0*sigma^2))
          *(P+m)=exp(razn/ksig); // numerator of Pmn equation
          sp+=*(P+m);// accumulate the denominator
      }
      
      // sp = sum{e^(|Xn - Ym|^2 / (-2.0*sigma^2))} + (2*pi*sigma^2)^(D/2)*w/(1-w)*M/N
      sp+=outlier_tmp; // add the outlier term to the denominator
      
      // Pt1[n] = 1 - ((2*pi*sigma^2)^(D/2)*w/(1-w)*M/N) / (sum{e^(|Xn - Ym|^2 / (-2.0*sigma^2))} + (2*pi*sigma^2)^(D/2)*w/(1-w)*M/N)
      *(Pt1+n)=1-outlier_tmp/ sp;
      
      // set temp x to current target point divided by (sum{e^(|Xn - Ym|^2 / (-2.0*sigma^2))} + (2*pi*sigma^2)^(D/2)*w/(1-w)*M/N)
      for (d=0; d < D; d++) {
       *(temp_x+d)=*(x+n+d*N)/ sp;
      }
      
      // for each point in the input set
      for (m=0; m < M; m++) {
         
          // P1[m] =  P[m] / ((sum{e^(|Xn - Ym|^2 / (-2.0*sigma^2))} + (2*pi*sigma^2)^(D/2)*w/(1-w)*M/N))
          *(P1+m)+=*(P+m)/ sp; // normalize the probability
          
          // for each dimension
          for (d=0; d < D; d++) {
          // Px[m][d] += temp_x[d]*P[m]
          *(Px+m+d*M)+= *(temp_x+d)**(P+m);
          }
          
      }
      
   *E +=  -log(sp);     
  }
  *E +=D*N*log(*sigma2)/2;
    
  
  free((void*)P);
  free((void*)temp_x);

  return;
}

/* Input arguments */
#define IN_x		prhs[0]
#define IN_y		prhs[1]
#define IN_sigma2	prhs[2]
#define IN_outlier	prhs[3]


/* Output arguments */
#define OUT_P1		plhs[0]
#define OUT_Pt1		plhs[1]
#define OUT_Px		plhs[2]
#define OUT_E		plhs[3]


/* Gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *sigma2, *outlier, *P1, *Pt1, *Px, *E;
  int     N, M, D;
  
  /* Get the sizes of each input argument */
  N = mxGetM(IN_x);
  M = mxGetM(IN_y);
  D = mxGetN(IN_x);
  
  /* Create the new arrays and set the output pointers to them */
  OUT_P1     = mxCreateDoubleMatrix(M, 1, mxREAL);
  OUT_Pt1    = mxCreateDoubleMatrix(N, 1, mxREAL);
  OUT_Px    = mxCreateDoubleMatrix(M, D, mxREAL);
  OUT_E       = mxCreateDoubleMatrix(1, 1, mxREAL);

    /* Assign pointers to the input arguments */
  x      = mxGetPr(IN_x);
  y       = mxGetPr(IN_y);
  sigma2       = mxGetPr(IN_sigma2);
  outlier    = mxGetPr(IN_outlier);

 
  
  /* Assign pointers to the output arguments */
  P1      = mxGetPr(OUT_P1);
  Pt1      = mxGetPr(OUT_Pt1);
  Px      = mxGetPr(OUT_Px);
  E     = mxGetPr(OUT_E);
   
  /* Do the actual computations in a subroutine */
  cpd_comp(x, y, sigma2, outlier, P1, Pt1, Px, E, N, M, D);
  
  return;
}


