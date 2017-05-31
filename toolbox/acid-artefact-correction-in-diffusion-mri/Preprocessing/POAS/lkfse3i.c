/* lkfse3i.c*/
/* 18/12/2012 */

/* [indtmp,wghttmp,ni] = lkfse3i(h(i),kappa(i),i,sdist,ng,vext,ni)*/

#include "mex.h"
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>

int maxSM(double m1, double m2) {
        if(m1>m2) 
            return (int)floor(m1); 
            else 
            return (int)floor(m2);
        }



void lkfse3i_main (double h, double kappa, double i4in, double *sdist, int ngin, double *vext, int niin,double *indtmp,double *wghttmp, double *niout)    
{
    int ni, i4, ng, ih1, ih2, ih3, i, j1, j2, j3, j4, sdistinx, indinx;
    float x1, x2, x3, k4, z, z1, vd2, vd3, h2; /* KT (26.04.2013): removed kap2, h2 */
    /*float *indtmp, *wghttmp;*/
/*    DD           = matrix( 1, N, 1, N);  
    eigenvectors = matrix( 1, N, 1, N);  
    eigenvalues  = vector( 1, N);  
    in_vox=0;
    out_vox=0;*/
    
    /*
    printf("h %f kapp: %f i4in %f ngin %i niin %i \n", h, kappa, i4in, ngin, niin); */   
    /*ni  = nii[0];*/
    vd2 = vext[0];
    vd3 = vext[1];
    h2  = h*h;
    ih1 = maxSM(1.0,h);/* JP need only h not 5*h*/
    ih2 = maxSM(1.0,h/vd2);
    ih3 = maxSM(1.0,h/vd3);
    i = 0;
    i4 = (int)i4in;
    ng = ngin;
    ni = niin;
    for (j4 = 0; j4 < ng; j4++){
        sdistinx = (i4) + j4*ng;
        k4  = sdist[sdistinx];
        z   = k4/kappa; /* KT (26.04.2013): need kappa and k4 instead of kap2 and k4*k4 here */
        if(z<=h){ /* KT (30.04.2013): need h here, not h2*/
            for (j1 = 0; j1 <= ih1; j1++){
                x1 = (float)j1;
                for (j2 = -ih2; j2 <= ih2; j2++){
                    x2 = (vd2*(float)j2);
                    for (j3 = -ih3; j3 <= ih3; j3++){
                        x3 = (vd3*(float)j3);
                        z1 = z+sqrt(x1*x1+x2*x2+x3*x3); /* KT (26.04.2013): included sqrt here*/        
                        if(z1<=h){ /* KT (30.04.2013): need h here, not h2*/
                            wghttmp[i]=1-z1*z1/h2; /* J.P. we use Epanechnicov which is k(x)=1-x^2 */
                            indinx = 0+i*5;
                            indtmp[indinx]=j1;
                            indinx = 1+i*5;
                            indtmp[indinx]=j2;                            
                            indinx = 2+i*5;
                            indtmp[indinx]=j3;                            
                            indinx = 3+i*5;
                            indtmp[indinx]=i4;                            
                            indinx = 4+i*5;
                            indtmp[indinx]=j4; /* was +1 KT+JP: 27.06.2013 */
                            i = i+1;
                        }/* end if z1 <= h2*/
                    }/* end j3*/
                }/* end j2*/
            }/* end j1*/
        }/* end if z <= h2*/
    }/* end j4*/
    niout[0] = (float)i;
  /*   printf("ni: %6f\n", niout[0]); */                        
    /*    niout[0] = 10.0;(double)i;, double *niout*/
}/* end main*/




void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* const int *dims;*/
    int out_ind[2], out_wght[2], out_ni[2], ni, ng;

    double h, kappa, *sdist, *vext, *indtmp, *wghttmp, i, *niout; /*, *niout*/ 
  

    if (nlhs > 4) mexErrMsgTxt ("Too many output arg."); /* to be adjusted */
    if (nrhs != 7) mexErrMsgTxt ("Wrong number of input arg.");  /* to be adjusted */
    
    
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
        mexErrMsgTxt ("Wrong input data type.");
    
    /*ndims = mxGetNumberOfDimensions (prhs[0]);*/
    
    /*if (ndims != 4) mexErrMsgTxt ("4d volume array req.");*/
    
        
    h       = mxGetScalar (prhs[0]); /* bandwidth */
    kappa	= mxGetScalar (prhs[1]); /* kappa */
    i		= mxGetScalar (prhs[2]); /* gradient number */
    sdist	= mxGetPr (prhs[3]); /* distant matrix for gradients */
	ng      = mxGetScalar (prhs[4]); /* number of gradients */
	vext	= mxGetPr (prhs[5]); /* vector of the length 2, relative voxel extensions */
	ni		= mxGetScalar (prhs[6]); /* maximal number of points */
	
   
       /* dims = mxGetDimensions (prhs[0]);8*/
    out_ind[0]     = 5; /* 1st dimension of ind*/
    out_ind[1]     = ni;    /*2nd dimension of ind */
    
    out_wght[0]     = 1;    /* 1st dimension of ni */
    out_wght[1]     = ni;   /*maxSMint(niout[0],ni);     1st dimension of ni */
    
    out_ni[0]      = 1;
    out_ni[1]      = 1;
    
	
    plhs[0] = mxCreateNumericArray (2, out_ind, mxDOUBLE_CLASS, mxREAL); 
    indtmp  = mxGetPr (plhs[0]);

    plhs[1] = mxCreateNumericArray (2, out_wght, mxDOUBLE_CLASS, mxREAL); 
    wghttmp = mxGetPr (plhs[1]);
    
    /*plhs[2] = mxCreateDoubleScalar(niout);*/
    
    plhs[2] = mxCreateNumericArray (2, out_ni, mxDOUBLE_CLASS, mxREAL); 
    niout   = mxGetPr (plhs[2]);/**/
    
    lkfse3i_main (h, kappa, i, sdist, ng, vext, ni,indtmp,wghttmp,niout);   
    
}
