/* ghfse3i.c */
/* used in gethseqfullse3 (getsphwghts.m)*/
/* 29/08/2013 */

/* [h(i,:),n] = ghfse3i(i,kstar,k456,ngrad,kappa,vext)*/

#include "mex.h"
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>

int maxInt(double m1, double m2) {
        if(m1<m2) 
            return m2; 
            else 
            return m1;
        }

void lkfse3i0(double h, double kappa, int i4, double *k456, int ng, double *vext, double *vred, int *nn){
/* we have dist=4 */ 
   int ih1,ih2,ih3,j1,j2,j3,j4,n;
   double x1,x2,x3,k4,z,z1,sw,sw2,wght,anz,h2,vd2,vd3;
   ih1 = (maxInt(1.e0,5.0e0*h));
   ih2 = (maxInt(1.e0,5.0e0*h/vext[0]));
   ih3 = (maxInt(1.e0,5.0e0*h/vext[1]));
   sw=0.e0;
   sw2=0.e0;
   h2 = h*h;
   vd2 = vext[0]*vext[0];
   vd3 = vext[1]*vext[1];
   n = 0;
   for (j4=0; j4<ng; j4++){ 
      k4 = k456[i4+ng*j4];
      z = k4/kappa;
      if(z<=h){
/* k4 not to large */
         for (j1=0; j1<=ih1; j1++){
/*   if j1>0  (-j1,-j2,-j3) gets the same weight, so count it twice */
            if(j1==0){
               anz=1.e0;
	    } else {
               anz=2.e0;
	    }
            x1 = j1;
            x1 = x1*x1;
            for (j2= -ih2; j2<=ih2; j2++){
               x2 = j2;
               x2 = x1 + vd2*x2*x2;
               if(x2 <= h2){
                  for (j3= -ih3; j3<=ih3; j3++){
                     x3 = j3;
                     z1 = x2+vd3*x3*x3;
                     z1=z+sqrt(z1);
                     if(z1<h){
                        wght= (1.e0-z1*z1/h2);
                        sw=sw+anz*wght;
                        wght=wght*wght;
                        sw2=sw2+anz*wght;
                        n=n+1;
		     }
	          }
	       }
	    }
	 }
      }
   }
/*   printf("sw %f sw2 %f \n", sw, sw2);*/
   *vred = sw*sw/sw2;
/*   printf("*vred %f  \n", *vred);*/
   *nn = n;
/*   printf(" *nn %i \n", *nn);*/
}

void ghfse3i_main(int i4, int kstar, double *k456, int ng, double kappa, double *vext, double *h, double *n){
/* compute bandwidth sequence for given kappa and gradients  
   lkfse3i0 computes weighting schemes and corresponding variance reduction
   k456,ng (input) contain auxiliary statistics(gradients) */
   int k,n0,maxn,nn;
   double hakt,hakt0,khakt,vr,ch,chk,vred,v0r;
   ch=1.25e0;
   hakt=1.e0;
   khakt= kappa/hakt;
/*   printf("hakt %f khakt %f i4 %d ng %d\n", hakt, khakt, i4, ng);*/
   lkfse3i0(hakt,khakt,i4,k456,ng,vext,&vr,&nn);
/*   printf("hakt %f vr %f \n", hakt, vr);*/
   chk=ch* vr;
   maxn = 1;
   for (k=0; k<kstar; k++){
      khakt= kappa/hakt;
      lkfse3i0(hakt,khakt,i4,k456,ng,vext,&vr,&nn);
/*  search for new hakt */  
      vred=vr/chk;
/*   printf("1: k %i hakt %f vr %f \n", k, hakt, vred);*/
      while (vred<1){
         hakt=hakt*1.05;
         khakt= kappa/hakt;
         lkfse3i0(hakt,khakt,i4,k456,ng,vext,&vr,&nn);
         vred=vr/chk;
/*   printf("2: hakt %f vr %f \n", hakt, vred);*/
      }
      while (vred>1){
         hakt0=hakt;
         v0r=vr;
         n0=nn;
         hakt=hakt/1.005;
         khakt= kappa/hakt;
         lkfse3i0(hakt,khakt,i4,k456,ng,vext,&vr,&nn);
         vred=vr/chk;
   /*printf("3: hakt %f vr %f \n", hakt, vred);*/
         if (vred<1) {
            hakt=hakt0;
            vr=v0r;
            nn=n0;
	 }
      }
      h[k] = hakt;
      chk=chk*ch;
      maxn=maxInt(maxn,nn);
      if(k==kstar){
      khakt= kappa/hakt;
      lkfse3i0(h[k],khakt,i4,k456,ng,vext,&vr,&nn);
   /*printf("4: hakt %f vr %f \n", hakt, vr);*/
     }
/*  number of positive weights for last step in n*/
   }
    *n = (float) maxn;
}


/* [h(i,:),n] = ghfse3i(i,kstar,k456,ngrad,kappa,vext)*/

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* const int *dims;*/
    int out_h[2], out_ni[2], i4, kstar, ng;
    double *k456, kappa, *vext, *h, *niout;
    if (nlhs > 2) mexErrMsgTxt ("Too many output arg."); /* to be adjusted */
    if (nrhs != 6) mexErrMsgTxt ("Wrong number of input arg.");  /* to be adjusted */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
        mexErrMsgTxt ("Wrong input data type.");
    i4 = mxGetScalar (prhs[0]);
    kstar = mxGetScalar (prhs[1]);
    k456 = mxGetPr (prhs[2]);
    ng = mxGetScalar (prhs[3]);
    kappa = mxGetScalar (prhs[4]);
    vext = mxGetPr (prhs[5]);
/*    printf( i4 %d kstar %d ng %d kappa %f\n", i4, kstar, ng, kappa);*/
   out_h[0] = 1;
    out_h[1] = kstar;
    plhs[0] = mxCreateNumericArray (2, out_h, mxDOUBLE_CLASS, mxREAL); 
    h  = mxGetPr (plhs[0]);

    out_ni[0]      = 1;
    out_ni[1]      = 1;
    plhs[1] = mxCreateNumericArray (2, out_ni, mxDOUBLE_CLASS, mxREAL); 
    niout   = mxGetPr (plhs[1]);/**/

    ghfse3i_main(i4, kstar, k456, ng, kappa, vext, h, niout);
}

