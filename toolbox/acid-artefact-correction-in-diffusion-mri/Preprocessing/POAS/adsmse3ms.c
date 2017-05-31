/* adsmse3ms.c*/
/* 28/8/2013 */

/* [nni, th, nni0, th0] = adsmse3ms(sb,s0,msth,msni,msth0,msni0,vmsth,vmsth0,mask,lambda,ws0,parind,parw,parnind,parind0,parw0,parnind0);
*/
#include "mex.h"
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>
#include "nrutil.h"
#include "nrutil.c"


void adsmse3ms(float *y, float *y0, float *th, float *ni, float *th0, float *ni0, float *fsi2, float *fsi02, float *mask, int ns, int n1, int n2, int n3, int ng, double lambda, double ws0, double *ind, double *w, int n, double *ind0, double *w0, int n0, float *thn, float *nin, float *th0n, float *ni0n){
/*
C  Multi-shell version (differs in dimension of th 
C  KL-distance based on all spheres and Gauss-approximation only
C  see ~polzehl/latex/1211_simplemetric/Rcode/Figuregaussapprox.r 
C  for approximative KL-distance between non-central chi distributions
C
C   perform adaptive smoothing on SE(3) multishell including s0
C   y  -  si images
C   y0 -  mean s0 image
C   th -  estimated/interpolated \E si on all shells (including 0 shell)
C   ni -  corresponding sum of weights
C   fsi2 - sd of noncentral chi with expectation th
C   th0 -  estimated/interpolated \E s0 and mean_g(\E si) on all other shells  
C   ni0 -  corresponding sum of weights
C   fsi02 - sd of noncentral chi with expectation th0
C   mask - head mask
C   ns   - number of shells (including 0 shell)
C   n1,n2,n3,ng - dimensions, number of gradients (bv!=0)
C   lambda - skale parameter
C   ws0  - relative weight for information from s0 images (should be in [0,1])
C   ind    - index vectors for si weighting schemes !! values in column 4 and 5 need to be in 0:(ng-1) !!
C   w    - corresponding weights
C   n    - number of si weights
C   ind0    - index vectors for s0 weighting schemes 
C   w0    - corresponding weights
C   n0    - number of s0 weights
C   thn   - new estimates of \E si
C   thn0   - new estimates of \E s0
C   nin   - new sum of weight for si
C   ni0n   - new sum of weight for s0
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w[i] for si images
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
*/
   int iind,n123,jind,i,i1,i2,i3,i4,iind4,j1,j2,j3,j4,jind4,j4n,k,cycle;
   double sz,z,sw0,swy0;
   double *sw=NULL, *swy=NULL, *fsi2i=NULL, *thi=NULL, *nii=NULL;

/*  First si - images */
   n123 = n1*n2*n3;
   sw      = dvector( 0, ng-1);  
   swy     = dvector( 0, ng-1);
   fsi2i   = dvector( 0, ns-1);
   thi   = dvector( 0, ns-1);
   nii   = dvector( 0, ns-1);
   for (iind=0; iind < n123; iind++){/* loop over voxel in R^3 */
      i1=iind%n1;
      i2=((iind-i1)/n1)%n2;
      i3=(((iind-i1)/n1-i2)/n2)%n3;/* i1,i2,i3 are the voxel coordinates */
      if(mask[iind]>0){ /* voxel is in the head mask */
	 for (i4=0; i4<ng; i4++){ 
            sw[i4]=0.e0;
            swy[i4]=0.e0;
	 }
         i4=-1;
         for (i=0; i<n; i++){/* loop over local neighborhood in SE(3) */
            cycle = 0;
            if(ind[3+5*i]!=i4){/* check if we are still considering gradient i4
	      if not assign correct i4 */
/*   by construction ind(4,.) should have same values consequtively */
               i4 = ind[3+5*i];
	       iind4 = iind*ns+i4*n123*ns;
	       for (k=0; k<ns; k++){/* get sd, theta and ni for voxel iind and gradient i4 */
                  fsi2i[k]=fsi2[k+iind4];
                  thi[k] = th[k+iind4];
                  nii[k] = ni[k+iind4]/lambda;
	       }
               nii[0]=ws0*nii[0];/* downweight s0 image */ 
	    }
            j1=i1+ind[5*i]; /* get coordinates of voxel j from description in ind */
            if((j1<0)|(j1>=n1)) cycle=1;
            j2=i2+ind[1+5*i];
            if((j2<0)|(j2>=n2)) cycle=1;
            j3=i3+ind[2+5*i];
            if((j3<0)|(j3>=n3)) cycle=1;
	    if(cycle==0) jind=j1+j2*n1+j3*n1*n2;
            if(cycle==0&mask[jind]>0){ /* if voxel in mask (and in the cube) 
	      compute adaptive weights and aggregate sum of weigths and weigted sum of observations in sw and swy */ 
               j4=ind[4+5*i]; 
	       j4n = jind+j4*n123;
	       jind4 = j4n*ns;
               sz=0.e0;
	       for (k=0; k<ns; k++){
                  z=(thi[k]-th[k+jind4]);
                  sz+=nii[k]*z*z/(fsi2[k+jind4]+fsi2i[k]);
	       }
               if(sz<1.e0){ 
                  z=w[i];
		  if(sz>.5) z=z*(2.e0-2.e0*sz);
                  sw[i4]=sw[i4]+z;
                  swy[i4]=swy[i4]+z*y[j4n];
	       }
	    }
	 }
/*  now the same for the opposite directions */
         for (i=0; i<n; i++){
            cycle = 0;
            if(ind[5*i]>0){
               if(ind[3+5*i]!=i4){
/*   by construction ind(4,.) should have same values consequtively */
                  i4 = ind[3+5*i];
	          iind4 = iind*ns+i4*n123*ns;
	          for (k=0; k<ns; k++){
                     fsi2i[k]=fsi2[k+iind4];
                     thi[k] = th[k+iind4];
                     nii[k] = ni[k+iind4]/lambda;
	          }
                  nii[0]=ws0*nii[0];
/*  first component corresponds to S0 image, ws0 is used to downweight its influence
C  when smoothing diffusion weighted data */
	       }
/*
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
*/
               j1=i1-ind[5*i];
               if(j1<0|j1>=n1) cycle=1;
               j2=i2-ind[1+5*i];
               if(j2<0|j2>=n2) cycle=1;
               j3=i3-ind[2+5*i];
               if(j3<0|j3>=n3) cycle=1;
	       if(cycle==0) jind=j1+j2*n1+j3*n1*n2;
               if(cycle==0&mask[jind]>0){
                  j4=ind[4+5*i];
   	          j4n = jind+j4*n123;
	          jind4 = j4n*ns;
                  sz=0.e0;
	          for (k=0; k<ns; k++){
                     z=(thi[k]-th[k+jind4]);
                     sz+=nii[k]*z*z/(fsi2[k+jind4]+fsi2i[k]);
	          }
                  if(sz<1.e0){ 
		     z=w[i];
		     if(sz>.5) z=z*(2.e0-2.e0*sz);
                     sw[i4]=sw[i4]+z;
                     swy[i4]=swy[i4]+z*y[j4n];
	          }
	       }
	    } 
         } 
/*  we now have everyting needed to compute the estimates for voxel i and all
  gradients (except s0)*/
         for (i4=0; i4<ng; i4++){
            thn[iind+i4*n123] = swy[i4]/sw[i4];
            nin[iind+i4*n123] = sw[i4];
         }
/*    now the s0 image in iind */
         sw0=0.e0;
         swy0=0.e0;
	 for (k=0; k<ns; k++){/* get sd, theta and ni for voxel iind  */
            thi[k] = th0[k+iind*ns];
            nii[k] = ni0[k+iind*ns]/lambda;
            fsi2i[k]=fsi02[k+iind*ns];
	 }
         for (i=0; i<n0; i++){/* loop over local neighborhood in R^3  */
            cycle = 0;
            j1=i1+ind0[3*i];
            if(j1<0|j1>=n1) cycle=1;
            j2=i2+ind0[1+3*i];
            if(j2<0|j2>=n2) cycle=1;
            j3=i3+ind0[2+3*i];
            if(j3<0|j3>=n3) cycle=1;
	    if(cycle==0) jind=j1+j2*n1+j3*n1*n2;
           if(cycle==0&mask[jind]>0){/* if voxel in mask (and in the cube) 
	      compute adaptive weights and aggregate sum of weigths and weigted sum of observations in sw and swy */
               sz=0.e0;
	       for (k=0; k<ns; k++){
                  z=(thi[k]-th0[k+jind*ns]);
                  sz+=nii[k]*z*z/(fsi02[k+jind*ns]+fsi2i[k]);
	       }
               if(sz<1.e0){
                  z=w0[i];
		  if(sz>.5) z=z*(2.e0-2.e0*sz);
                  sw0+=z;
                  swy0+=z*y0[jind];
	       }
	    }
	 }
/*  now opposite directions */
         for (i=0; i<n0; i++){
           cycle = 0;
           if(ind0[3*i]>0){
/*
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
*/
               j1=i1-ind0[3*i];
               if(j1<0|j1>=n1) cycle=1;
               j2=i2-ind0[1+3*i];
               if(j2<0|j2>=n2) cycle=1;
               j3=i3-ind0[2+3*i];
               if(j3<0|j3>=n3) cycle=1;
	       if(cycle==0) jind=j1+j2*n1+j3*n1*n2;
               if(cycle==0&mask[jind]>0){
                  sz=0.e0;
	          for (k=0; k<ns; k++){
                     z=(thi[k]-th0[k+jind*ns]);
                     sz+=nii[k]*z*z/(fsi02[k+jind*ns]+fsi2i[k]);
	          }
                  if(sz<=1.e0){
                     z=w0[i];
   	   	     if(sz>.5) z=z*(2.e0-2.e0*sz);
                     sw0+=z;
                     swy0+=z*y0[jind];
	          }
	       }
	    }
         }
/*  we now have everyting needed to compute the estimates for voxel i and zero gradient (S0)*/
         th0n[iind] = swy0/sw0;
         ni0n[iind] = sw0;
      }
   }
   free_dvector(sw, 0, ng-1);  
   free_dvector(swy, 0, ng-1);
   free_dvector(fsi2i, 0, ns-1);
   free_dvector(thi, 0, ns-1);
   free_dvector(nii, 0, ns-1);
}



/* [nni, th, nni0, th0] = adsmse3ms(sb,s0,msth,msni,msth0,msni0,vmsth,vmsth0,mask,lambda,ws0,parind,parw,parnind, parind0,parw0,parnind0);
*/


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int ndims;
    int i, ns, ng, n1, n2, n3, parnind, parnind0;
    const int *dims;
    int out_th[4], out_th0[3];
    double lambda, ws0;

    float *si=NULL, *s0=NULL, *th=NULL, *th0=NULL, *ni=NULL, *ni0=NULL, *vth=NULL, *vth0=NULL, *mask=NULL;
    double *parind=NULL, *parw=NULL, *parind0=NULL, *parw0=NULL;
    float *thout=NULL, *niout=NULL, *th0out=NULL, *ni0out=NULL;
  

    if (nlhs > 4) mexErrMsgTxt ("Too many output arg."); /* to be adjusted */
    if (nrhs != 17) mexErrMsgTxt ("Wrong number of input arg.");  /* to be adjusted */
    
    
    dims        = mxGetDimensions (prhs[2]);    
    ns          = dims[0];
    n1          = dims[1];
    n2          = dims[2];
    n3          = dims[3];
    ng          = dims[4];
    
    
    si          = (float*) mxGetPr (prhs[0]); /* original dw images*/
    s0          = (float*) mxGetPr (prhs[1]); /* original b0 images*/
    th          = (float*) mxGetPr (prhs[2]); /* estimated and interpolated dw images on all shells*/
    ni          = (float*) mxGetPr (prhs[3]); /* sum of weights for dw images*/
    th0         = (float*) mxGetPr (prhs[4]); /* estimated and interpolated b0 images on all shells*/
    ni0         = (float*) mxGetPr (prhs[5]); /* sum of weights for b0 images*/
    vth         = (float*) mxGetPr (prhs[6]); /* variance correction factors*/
    vth0        = (float*) mxGetPr (prhs[7]); /* variance correction factors*/
    mask        = (float*) mxGetPr (prhs[8]); /* brain mask */
    lambda      = mxGetScalar (prhs[9]); /* scale factor for statistical penalty */
    ws0         = mxGetScalar (prhs[10]); /* downweighting of s0 */
    parind      = mxGetPr (prhs[11]); /* index field - 5xn. - integer!!!*/
    parw        = mxGetPr (prhs[12]); /* weights n elements*/
    parnind     = (int) mxGetScalar (prhs[13]); /* integer scalar: "n" */
    parind0      = mxGetPr (prhs[14]); /* index field - 3xn. - integer!!!*/
    parw0        = mxGetPr (prhs[15]); /* weights n elements*/
    parnind0     = (int) mxGetScalar (prhs[16]); /* integer scalar: "n" */

	
    
       /* dims = mxGetDimensions (prhs[0]);8*/
    out_th[0]   = n1; /* 1st dimension of ind*/
    out_th[1]   = n2;   /* 2nd dimension of ind */
    out_th[2]   = n3; /* 1st dimension of ind*/
    out_th[3]   = ng;   /* 2nd dimension of ind */
    
    plhs[0]     = mxCreateNumericArray (4, out_th, mxSINGLE_CLASS, mxREAL); 
    niout       = (float*) mxGetPr (plhs[0]);

    plhs[1]     = mxCreateNumericArray (4, out_th, mxSINGLE_CLASS, mxREAL); 
    thout       = (float*) mxGetPr (plhs[1]);
    
       /* dims = mxGetDimensions (prhs[0]);8*/
    out_th0[0]   = n1; /* 1st dimension of ind*/
    out_th0[1]   = n2;   /* 2nd dimension of ind */
    out_th0[2]   = n3; /* 1st dimension of ind*/
    
    plhs[2]     = mxCreateNumericArray (3, out_th0, mxSINGLE_CLASS, mxREAL); 
    ni0out       = (float*) mxGetPr (plhs[2]);

    plhs[3]     = mxCreateNumericArray (3, out_th0, mxSINGLE_CLASS, mxREAL); 
    th0out       = (float*) mxGetPr (plhs[3]);

    adsmse3ms(si, s0, th, ni, th0, ni0, vth, vth0, mask, ns, n1, n2, n3, ng, lambda, ws0, parind, parw, parnind, parind0, parw0, parnind0, thout, niout, th0out, ni0out);    
}
