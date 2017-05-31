/* DTI2EVEW.c */

/* [FA, EV1, EV2, EV3, EW, outnr] = DTI2EVEW (VOL) */

#include "mex.h"
#include <stdlib.h> 
#include <math.h>
/*#include <stddef.h>*/
#include <stdio.h>
#include "nrutil.h"
#include "nrutil.c"
#include "jacobi.c"
#include "gaussj_mod.c"

#define N 3
#define ZEILEN 6
#define LIMIT_ADC 0.005
#define THRESHOLD 5
#define CONS 1.1 
struct point
{
    int xa;
    int n;
};
 
void GVtoTens (double **TV, double *GV, int NG)
 {
    int j, GVpos_X, GVpos_Y, GVpos_Z;
    double x, y, z, n;
    
    for (j = 0; j < NG;j++)
    {   
        x            = 0.0;
        y            = 0.0;
        z            = 0.0;
        n            = 0.0;
        GVpos_X      = j * N;
        GVpos_Y      = 1 + j * N;
        GVpos_Z      = 2 + j * N;
        x            = GV[GVpos_X];
        y            = GV[GVpos_Y];
        z            = GV[GVpos_Z];
        n            = sqrt(x*x+y*y+z*z);
        TV[(j+1)][1] = x*x/n;
        TV[(j+1)][2] = y*y/n;
        TV[(j+1)][3] = z*z/n;
        TV[(j+1)][4] = 2*x*y/n;
        TV[(j+1)][5] = 2*x*z/n;
        TV[(j+1)][6] = 2*y*z/n;
    }
 }

 void TVtoXM (double **XM, double **TV, double *SIG, int NG) 
 {
    int i, j, k;
    double sum;
    
    for (i=1; i <= ZEILEN; i++)
    {
        for (j=1; j <= ZEILEN; j++)
        {
            sum = 0.0;
            for (k=1; k <= NG; k++)
            {
                sum = TV[k][i] * TV[k][j] * SIG[k] + sum;
            }
            XM[i][j] = sum;
   
        }
    }
 }

 void NteSPALTE_XM (double **VEC, double **XM, double N_XM) 
 {
    int i, NX;
    NX = N_XM;
    for (i=1; i <= ZEILEN; i++)
    {
        VEC[i][1] = XM[i][NX];
    }
 }


int max1(float *egv)
    {  
        if (egv[1] > egv[2])
        {
            if (egv[1] > egv[3]) return 1; 
            else return 3;
        } 
		else 
		{
		    if (egv[2] > egv[3]) return 2;
            else return 3;
        }
	}
		       
int max2(float *egv)
    {  
        if (egv[1] > egv[2])
        {
            if (egv[1] > egv[3]) 
            {
                if (egv[2] > egv[3])  return 2; 
                else return 3;
            }
			else return 1;
        } 
		else 
		{
		    if (egv[2] > egv[3]) 
		    {
		        if (egv[1] > egv[3])  return 1; 
			    else return 3;
			}
            else return 2;
        }
   }
    
int max3(int m1, int m2) {
        if((m1-m2)*(m1-m2)==1)
        {  
            if((m1==3) || (m2==3)) return 1; 
            else return 3;
        }
        else return 2;
   }

void VoxelPos ( double *FA, double *iVOL, const int *dims, double *EV1, double *EV2, double *EV3, double *EW, double *outlV)
{
    long int pos_iVOL , pos, pos_E, out_vox; /*relation between position in array (DT/iVOL) and voxel-space */
    
    struct point i;
    int j, in_vox, num_rot,  eig1, eig2, eig3;
    float av, dM[6];
    float **eigenvectors, *eigenvalues, **DD;
    DD           = matrix( 1, N, 1, N);  /* must be a symmetric real matrix */
    eigenvectors = matrix( 1, N, 1, N);  /* eigenvectors as columns; first index classifies the range of the eigenvector*/
    eigenvalues  = vector( 1, N);  
    in_vox=0;
    out_vox=0;
    
    for (i.xa = 0; i.xa < dims[0]; i.xa++)
    {
        pos         = i.xa;   /* x*y*z-Position */    
        outlV[pos]  = 0;

        for (i.n = 0; i.n < dims[1]; i.n++) 
        {   
            pos_iVOL = i.xa + i.n*dims[0]; 
            dM[i.n]  = iVOL[pos_iVOL];
        }

        if((dM[0]+dM[1]+dM[2]) >0) 
        {
            DD[1][1]   = dM[0];
            DD[2][2]   = dM[1];
            DD[3][3]   = dM[2];
            DD[1][2]   = DD[2][1] = dM[3];
            DD[1][3]   = DD[3][1] = dM[4];
            DD[2][3]   = DD[3][2] = dM[5];

            jacobi(DD, N, eigenvalues, eigenvectors, &num_rot);
            av          = (DD[1][1] + DD[2][2] + DD[3][3])/3;
            FA[pos]     = 1000.0 * sqrt( 3 * ((eigenvalues[1]-av) * (eigenvalues[1]-av) + (eigenvalues[2]-av) * (eigenvalues[2]-av) + (eigenvalues[3]-av) * (eigenvalues[3]-av) )) / sqrt( 2.0 * ( eigenvalues[1] * eigenvalues[1] + eigenvalues[2] * eigenvalues[2] + eigenvalues[3] * eigenvalues[3]) );
            
            eig1        = max1(eigenvalues);
            eig2        = max2(eigenvalues);
            eig3        = max3(eig1,eig2);

            pos_E      = pos;
            EW[pos_E]  = eigenvalues[eig1];
            EV1[pos_E]  = eigenvectors[1][eig1];
            EV2[pos_E]  = eigenvectors[2][eig1];
            EV3[pos_E]  = eigenvectors[3][eig1];
            pos_E      = pos + 1*dims[0];
            EW[pos_E]  = eigenvalues[eig2];
            EV1[pos_E]  = eigenvectors[1][eig2];
            EV2[pos_E]  = eigenvectors[2][eig2];
            EV3[pos_E]  = eigenvectors[3][eig2];
            pos_E      = pos + 2*dims[0];
            EW[pos_E]  = eigenvalues[eig3];
            EV1[pos_E]  = eigenvectors[1][eig3];
            EV2[pos_E]  = eigenvectors[2][eig3];
            EV3[pos_E]  = eigenvectors[3][eig3];

            in_vox++ ;   
        }
        else
        {
            FA[pos]     = 0.0;
            
            pos_E       = pos;
            EW[pos_E]   = 0.0;
            EV1[pos_E]   = 0.0;
            EV2[pos_E]   = 0.0;
            EV3[pos_E]   = 0.0;
            pos_E       = pos + 1*dims[0];
            EW[pos_E]   = 0.0;
            EV1[pos_E]   = 0.0;
            EV2[pos_E]   = 0.0;
            EV3[pos_E]   = 0.0;
            pos_E       = pos + 2*dims[0];
            EW[pos_E]   = 0.0;
            EV1[pos_E]   = 0.0;
            EV2[pos_E]   = 0.0;
            EV3[pos_E]   = 0.0;


        }
        
    }/* i.x */
    printf("Number of voxels with outlier: %d\n", out_vox);
}




void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int ndims;
    const int *dims;
    int out_dims[2];
    int out3_dims[2];
    int out4_dims[2];

    double *FA, *VOL, *EV1, *EV2, *EV3, *EW, *outlV; 
  

    if (nlhs > 6) mexErrMsgTxt ("Too many output arg.");
    if (nrhs != 1) mexErrMsgTxt ("Wrong number of input arg.");
    
    
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
        mexErrMsgTxt ("Wrong input data type.");
    
    ndims = mxGetNumberOfDimensions (prhs[0]);
    
    if (ndims != 2) mexErrMsgTxt ("4d volume array req.");
    
    
    dims = mxGetDimensions (prhs[0]);
    out_dims[0] = dims[0];
    out_dims[1] = 1;
    out3_dims[0] = dims[0];
    out3_dims[1] = 3;
    out4_dims[0] = dims[0];
    out4_dims[1] = 6;
    
    VOL     = mxGetPr (prhs[0]);
   
    plhs[0] = mxCreateNumericArray (ndims, out_dims, mxDOUBLE_CLASS, mxREAL); /* residual error */
    FA      = mxGetPr (plhs[0]);

    plhs[1] = mxCreateNumericArray (ndims, out3_dims, mxDOUBLE_CLASS, mxREAL); /* Eigenvectors */
    EV1    = mxGetPr (plhs[1]);    
    plhs[2] = mxCreateNumericArray (ndims, out3_dims, mxDOUBLE_CLASS, mxREAL); /* Eigenvectors */
   EV2    = mxGetPr (plhs[2]);
    plhs[3] = mxCreateNumericArray (ndims, out3_dims, mxDOUBLE_CLASS, mxREAL); /* Eigenvectors */
    EV3    = mxGetPr (plhs[3]);
    plhs[4] = mxCreateNumericArray (ndims, out3_dims, mxDOUBLE_CLASS, mxREAL); /* Eigenvectors */
    EW    = mxGetPr (plhs[4]);
    plhs[5] = mxCreateNumericArray (ndims, out_dims, mxDOUBLE_CLASS, mxREAL); /* Map with outliers */
    outlV   = mxGetPr (plhs[5]);
    
    
    VoxelPos (FA,VOL, dims, EV1, EV2, EV3, EW,outlV);    
}
