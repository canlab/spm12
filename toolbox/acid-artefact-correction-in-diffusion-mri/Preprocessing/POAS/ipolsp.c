/* ipolsph used in interpolatesphere0 */
/* 18/12/2012 */

/*[mstheta(:,mask,:),msni(:,mask,:)] = ipolsp(theta(mask,:),th0(mask),ni(mask,:),ni0(mask),nmask,ng,n3g{4},n3g{5},nbv); */

#include "mex.h"
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>

void ipolsp_main (float *theta, float *th0, float *ni, float *ni0, int n, int ng, float *gind, float* gw, int nbv, float *msth, float *msni)
{
/*   interpolate values of theta on spheres where it was not observed */
/*      integer n,ng,nbv,nbvp1,gind(3,nbv,ng)
      real*8 theta(n,ng),th0(n),ni(n,ng),ni0(n),gw(3,nbv,ng),
     1       msth(nbvp1,n,ng),msni(nbvp1,n,ng) */
      int i,j,k,i1,i2,i3,ip1, nbvp1, indjk; 
      float w1,w2,w3;
      nbvp1 = nbv + 1;
      for (j=0;j<ng;j++){
	  for (k=0;k<n;k++){
            indjk = k*nbvp1+j*n*nbvp1;
            msth[indjk] = th0[k];
            msni[indjk] = ni0[k];
	    for (i=0;i<nbv;i++){
               ip1 = i+1;
               i1 = gind[i*3+j*nbv*3];
	       if(i1==j){
		  msth[ip1+indjk] = theta[k+j*n];
		  msni[ip1+indjk] = ni[k+j*n];
	       } else {
/*  different shell need to interpolate */
                  i2 = gind[1+i*3+j*nbv*3];
                  i3 = gind[2+i*3+j*nbv*3];
                  w1 = gw[i*3+j*nbv*3];
                  w2 = gw[1+i*3+j*nbv*3];
                  w3 = gw[2+i*3+j*nbv*3];
                  msth[ip1+indjk] = theta[k+i1*n]*w1+theta[k+i2*n]*w2+theta[k+i3*n]*w3;
                  msni[ip1+indjk] = 1.e0/(w1/ni[k+i1*n]+w2/ni[k+i2*n]+w3/ni[k+i3*n]);
	       }
	    }
	 }
      }
}      

/*[mstheta(:,mask,:),msni(:,mask,:)] = ipolsp(theta(mask,:),th0(mask),ni(mask,:),ni0(mask),nmask,ng,n3g{4},n3g{5},nbv); */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int out_msth[3], n, ng, nbv;
    float *theta, *th0, *ni, *ni0, *gw, *msth, *msni, *gind;
    if (nlhs > 2) mexErrMsgTxt ("Too many output arg."); /* to be adjusted */
    if (nrhs != 9) mexErrMsgTxt ("Wrong number of input arg.");  /* to be adjusted */
    theta = (float*) mxGetData (prhs[0]);
    th0 = (float*) mxGetData (prhs[1]);
    ni = (float*) mxGetData (prhs[2]);
    ni0 = (float*) mxGetData (prhs[3]);
    n = mxGetScalar (prhs[4]);
    ng = mxGetScalar (prhs[5]);
    gind = (float*) mxGetData (prhs[6]);
    gw = (float*) mxGetData (prhs[7]);
    nbv = mxGetScalar (prhs[8]);
    
    out_msth[0]     = nbv+1; /* 1st dimension of msth*/
    out_msth[1]     = n;    /*2nd dimension of msth */
    out_msth[2]     = ng;    /*2nd dimension of msth */
    
    plhs[0] = mxCreateNumericArray (3, out_msth, mxSINGLE_CLASS, mxREAL); 
    msth  = (float*) mxGetData (plhs[0]);
    plhs[1] = mxCreateNumericArray (3, out_msth, mxSINGLE_CLASS, mxREAL); 
    msni  = (float*) mxGetData (plhs[1]);

    ipolsp_main(theta, th0, ni, ni0, n, ng, gind, gw, nbv, msth, msni);
}