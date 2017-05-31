/* interpol.c */
/* used in POAS_v02.m*/
/* 29/08/2013 */

/*        vmsth = linterpol(th,table,in+1,out+1);
*/

#include "mex.h"
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>
#include "nrutil.h"
#include "nrutil.c"

void linterpol(float *x, double *table, int in, int out, int n, float *fx){
/* in should be 0 or 1 */
   int i,j;
   double xi, x1, x2, fx1, fx2;
   for ( i=0; i<n; i++){/* loop over voxel */
      xi = x[i];
      if(xi>table[in*122]){/* x larger than first table entry */
         j = 1;
         if(xi<table[121+in*122]){/* x smaller than last table entry */
            while (xi > table[j+in*122]){
	       j=j+1;
            }
            x1 = table[j-1+in*122];
            fx1 = table[j-1+out*122];
            x2 = table[j+in*122];
            fx2 = table[j+out*122];
	    fx[i] = fx1 + (fx2-fx1)/(x2-x1)*(xi-x1);
         } else {/* x >= than last table entry */
  	    if(out<2){/* we have ncp \approx mu */
	       fx[i]=xi;
	    } else {/* we have sd and var \approx 1 */
	       fx[i] = 1;
	    } 
         }	
      } else {/* x <= first table entry use minimal fx */
	 fx[i] = table[out*122];
      }
   }/* end loop over voxel */
}
/*        vmsth = linterpol(th,table,in+1,out+1);
*/

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* const int *dims;*/
    const int *dims;
    int out_fx[2], n, in, out, ldims, i;
    float *x, *fx;
    double *table;
    if (nlhs > 1) mexErrMsgTxt ("Too many output arg."); /* to be adjusted */
    if (nrhs != 4) mexErrMsgTxt ("Wrong number of input arg.");  /* to be adjusted */
    ldims       = mxGetNumberOfDimensions(prhs[0]);
    dims        = mxGetDimensions (prhs[0]);    
    n = 1;
    for(i=0; i<ldims; i++){ 
       n = n* dims[i];
    }/* n contains number of voxel x number of grads*/
    x = (float*) mxGetPr (prhs[0]);
    table = mxGetPr (prhs[1]);
    in = mxGetScalar (prhs[2])-1;
    out = mxGetScalar (prhs[3])-1;
    out_fx[0]     = 1;    /* 1st dimension of bg and bghat */
    out_fx[1]     = n;   /* 2nd dimension of bg and bghat */
 
    plhs[0] = mxCreateNumericArray (2, out_fx, mxSINGLE_CLASS, mxREAL);
    fx = (float*) mxGetPr (plhs[0]);
    linterpol(x, table, in, out, n, fx);
}
