/* lkfse3i.c*/
/* 18/12/2012 */

/* [parind0, parw0, parnind0] = lkfulls0(hakt0,vext,nind);*/

#include "mex.h"
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>


int maxInt(double m1, double m2) {
        if(m1>m2) 
            return m1; 
            else 
            return m2;
        }

void lkfulls0_main(double h, double *vext, int n, double *ind, double *wght,  double *nout){
   int ih1,ih2,ih3,i,j1,j2,j3;
   double h2,x1,x2,x3,z1,vd2,vd3;
   vd2 = vext[0];
   vd3 = vext[1];
   ih1 = maxInt(1.e0,h);
   ih2 = maxInt(1.e0,h/vd2);
   ih3 = maxInt(1.e0,h/vd3);
/*   printf("ih1 %d ih2 %d ih3 %d n %d \n",ih1,ih2,ih3,n);*/
   h2 = h*h;
   i = 0;
   for ( j1=0; j1<=ih1; j1++){
      x1 = j1;
      for ( j2=-ih2; j2<=ih2; j2++){
         x2 = vd2*j2;
         for ( j3=-ih3; j3<=ih3; j3++){
	    x3 = vd3*j3;
            z1 = x1*x1+x2*x2+x3*x3;
	    if(z1<h2){
	       if(i>=n){
		  printf("number of weights exceeds limit %d in lkfulls0\n",n);
		  i=n-1;
	       }
               wght[i]= (1.e0-z1/h2);
               ind[3*i] = j1;
               ind[3*i+1] = j2;
               ind[3*i+2] = j3;
               i = i+1;
	    }
	 }
      } 
   }
   *nout = (float)i;
}
/* end main*/




void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* const int *dims;*/
    int out_ind[2], out_wght[2], out_ni[2], ni;

    double h, *vext, *indtmp, *wghttmp,  *nout; /*, *nout*/ 
  

    if (nlhs > 3) mexErrMsgTxt ("Too many output arg."); /* to be adjusted */
    if (nrhs != 3) mexErrMsgTxt ("Wrong number of input arg.");  /* to be adjusted */
    
    
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
        mexErrMsgTxt ("Wrong input data type.");
    
    
        
    h       = mxGetScalar (prhs[0]); /* bandwidth */
    vext    = mxGetPr (prhs[1]); /* vector of the length 2, relative voxel extensions */
    ni      = (int) mxGetScalar (prhs[2]); /* maximal number of points */
	
/*    printf("h %f vext %f %f ni %d \n",h, vext[0], vext[1] ,ni);*/
  
       /* dims = mxGetDimensions (prhs[0]);8*/
    out_ind[0]     = 3; /* 1st dimension of ind*/
    out_ind[1]     = ni;    /*2nd dimension of ind */
    
    out_wght[0]     = 1;    /* 1st dimension of ni */
    out_wght[1]     = ni;   /*maxSMint(nout[0],ni);     1st dimension of ni */
    
    out_ni[0]      = 1;
    out_ni[1]      = 1;
    
	
    plhs[0] = mxCreateNumericArray (2, out_ind, mxDOUBLE_CLASS, mxREAL); 
    indtmp  = mxGetPr (plhs[0]);

    plhs[1] = mxCreateNumericArray (2, out_wght, mxDOUBLE_CLASS, mxREAL); 
    wghttmp = mxGetPr (plhs[1]);
    
    /*plhs[2] = mxCreateDoubleScalar(nout);*/
    
    plhs[2] = mxCreateNumericArray (2, out_ni, mxDOUBLE_CLASS, mxREAL); 
    nout   = mxGetPr (plhs[2]);/**/
    
    lkfulls0_main(h, vext, ni,indtmp,wghttmp,nout);   
    
}
