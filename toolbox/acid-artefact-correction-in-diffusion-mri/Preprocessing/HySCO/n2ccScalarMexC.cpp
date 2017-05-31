/*
 * (c) Lars Ruthotto 2013
 * http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
 *
 * CPP code for  n2ccScalar
 */

#include <math.h>
#include <mex.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

void nodalToCenterTwoD   (double *Pu,const double *u, const double *m);
void centerToNodalTwoD   (double *Pu,const double *u, const double *m);
void nodalToCenterThreeD (double *Pu,const double *u, const double *m);
void centerToNodalThreeD (double *Pu,const double *u, const double *m);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int dim, i;
    if (nrhs<2)
        mexErrMsgTxt("Number of arguments must be 2!");
    
    // Get input
    const double* u = static_cast<double*>(mxGetData(prhs[0]));
    const double* m = static_cast<double*>(mxGetData(prhs[1]));
    // get dimension, h and number of elements N
    dim = mxGetN(prhs[1]);
    double Ncc = 1, Nn=1;;
    
    for (int i=0; i < dim; i++)
    {  
        Ncc *= m[i];
        Nn *= m[i]+1;
    }
    
    const int N = mxGetM(prhs[0]);
    
    if (N==Ncc)
    {
        mwSize dims[2]; dims[0] = Nn; dims[1] = 1;
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        double* Pu = static_cast<double*>(mxGetData(plhs[0]));
        switch (dim) {
            case 2:
                centerToNodalTwoD(Pu,u,m);
                break;
            case 3:
                centerToNodalThreeD(Pu,u,m);
                break;
            default:
                break;
        }
    }
    else if(N==Nn)
    {
        mwSize dims[2]; dims[0] = Ncc; dims[1] = 1;
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        double* Pu = static_cast<double*>(mxGetData(plhs[0]));
         switch (dim) {
            case 2:
                nodalToCenterTwoD(Pu,u,m);
                break;
            case 3:
                nodalToCenterThreeD(Pu,u,m);
                break;
            default:
                break;
        }
   }
    else
    {
        mexErrMsgTxt("Unsupported input!");
    }

  }

void nodalToCenterTwoD (double *Pu,const double *u, const double *m)
{
    const int m1 = (int) m[0];
    const int m2 = (int) m[1];
    const int inc_x = 1;
    const int inc_y = m1;
    const int inc_yn = m1+1;
   
    const double *rptr;
    double *wptr;
    int i1,i2;
    #pragma omp parallel for default(shared) private(rptr,wptr,i1, i2)
    for (i2=0; i2<m2; i2++){
        for (i1=0; i1<m1; i1++){
                wptr = Pu + i1 + i2*(inc_y);
                rptr = u  + i1 + i2* inc_yn;
               // add A
                *wptr = *rptr;
                
                // add B
                rptr += inc_x;
                *wptr += *rptr;
                
                // add D
                rptr +=  (inc_yn);
                *wptr += *rptr;
                  
                // add C
                rptr  -= inc_x;
                *wptr += *rptr;
                
                // divide by four
                *wptr *= 0.25;
            }
        }
    
}
void centerToNodalTwoD (double *Pu,const double *u, const double *m)
{
    const int m1 = (int) m[0];
    const int m2 = (int) m[1];
    const int inc_x = 1;
    const int inc_y = m1;
    const int inc_yn = m1+1;
    
   
    double p[4];
 
    int i1,i2;
    for ( i2=0; i2<m2+1; i2++){
        for ( i1=0; i1<m1+1; i1++){
               // divide by four
                Pu[i1 + i2*inc_yn]  = ((i1>0) && (i2>0))             ? u[ (i1-1) + (i2-1)*inc_y ]  : 0.0 ; 
                Pu[i1 + i2*inc_yn] += ((i1>0) && (i2<(m2)) )       ? u[ (i1-1) + i2    *inc_y ]  : 0.0 ; 
                Pu[i1 + i2*inc_yn] += ((i2>0) && (i1<(m1)) )       ? u[ i1     + (i2-1)*inc_y ]  : 0.0 ; 
                Pu[i1 + i2*inc_yn] += ((i1<(m1))   && (i2<(m2))) ? u[ i1     + i2    *inc_y ]  : 0.0 ;
                
                 // add to A
                Pu[i1 + i2*inc_yn] *= 0.25;
             }
        }
   
}
void nodalToCenterThreeD (double *Pu,const double *u, const double *m)
{
    const int m1 = (int) m[0];
    const int m2 = (int) m[1];
    const int m3 = (int) m[2];
    const int inc_x = 1;
    const int inc_y = m1;
    const int inc_yn = m1+1;
    const int inc_z = m1*m2;
    const int inc_zn = (m1+1)* (m2+1);
 
    const double *rptr;
    double *wptr;
    int i1,i2,i3;
    #pragma omp parallel for default(shared) private(rptr,wptr,i1, i2,i3)
    for (i3=0; i3<m3; i3++){
        for (i2=0; i2<m2; i2++){
            for (i1=0; i1<m1; i1++){
                wptr = Pu + i1 + i2*(inc_y) + i3*(inc_z);
                rptr = u  + i1 + i2* inc_yn + i3*(inc_zn);

                // add A
                *wptr = *rptr;
                
                // add B
                rptr += inc_x;
                *wptr += *rptr;
                
                // add D
                rptr +=  inc_yn;
                *wptr += *rptr;
                
                // add C
                rptr  -= inc_x;
                *wptr += *rptr;
                
                // add G
                rptr  += inc_zn;
                *wptr += *rptr;
                
                // add H
                rptr += inc_x;
                *wptr += *rptr;
                
                // add F
                rptr  -= (inc_yn);
                *wptr += *rptr;
                
                // add E
                rptr  -= inc_x;
                *wptr += *rptr;
                // divide by eight
                
                *wptr *= 0.125;
            }
        }
    }

}
void centerToNodalThreeD (double *Pu,const double *u, const double *m)
{
    const int m1 = (int) m[0];
    const int m2 = (int) m[1];
    const int m3 = (int) m[2];
    const int inc_x = 1;
    const int inc_y = m1;
    const int inc_yn = m1+1;
    const int inc_z = m1*m2;
    const int inc_zn = (m1+1)*(m2+1);
    
    double *wptr;
    
    int i1,i2,i3;
    #pragma omp parallel for default(shared) private(wptr,i1, i2,i3)
    for (int i3=0; i3<m3+1; i3++){
        for (int i2=0; i2<m2+1; i2++){
            for (int i1=0; i1<m1+1; i1++){
                
                wptr = Pu + i1 + i2*inc_yn + i3*inc_zn;

                // divide by four
                *wptr += ((i1==0)  || (i2==0)  || (i3==0)) ? 0: u[ (i1-1) +(i2-1)*inc_y +(i3-1)*inc_z];
                *wptr += ((i1==m1) || (i2==0)  || (i3==0)) ? 0: u[ (i1)   +(i2-1)*inc_y +(i3-1)*inc_z];
                *wptr += ((i1==0)  || (i2==m2) || (i3==0)) ? 0: u[ (i1-1) +(i2)*inc_y   +(i3-1)*inc_z];
                *wptr += ((i1==m1) || (i2==m2) || (i3==0)) ? 0: u[ (i1)   +(i2)*inc_y   +(i3-1)*inc_z];
                
                *wptr += ((i1==0)  || (i2==0)  || (i3==m3)) ? 0: u[ (i1-1) +(i2-1)*inc_y +(i3)*inc_z];
                *wptr += ((i1==m1) || (i2==0)  || (i3==m3)) ? 0: u[ (i1)   +(i2-1)*inc_y +(i3)*inc_z];
                *wptr += ((i1==0)  || (i2==m2) || (i3==m3)) ? 0: u[ (i1-1) +(i2)*inc_y   +(i3)*inc_z];
                *wptr += ((i1==m1) || (i2==m2) || (i3==m3)) ? 0: u[ (i1)   +(i2)*inc_y   +(i3)*inc_z];
         
                *wptr *= 0.125;
       
            }
        }
    }
    
}

/*
    (c) Lars Ruthotto and Jan Modersitzki 2013

    This file is part of HySCO (Version 1.0, 2013/03/28)
                           -  Hyperelastic Susceptibility Artefact Correction for DTI

    
    HySCO is free but copyright software, distributed under the terms of the 
    GNU General Public Licence as published by the Free Software Foundation 
    (Version 3, 29 June 2007) http://www.gnu.org/licenses/gpl.html

 
    This code is provided "as is", without any warranty of any kind, either
    expressed or implied, including but not limited to, any implied warranty
    of merchantibility or fitness for any purpose. In no event will any party
    who distributed the code be liable for damages or for any claim(s) by
    any other party, including but not limited to, any lost profits, lost
    monies, lost data or data rendered inaccurate, losses sustained by
    third parties, or any other special, incidental or consequential damages
    arising out of the use or inability to use the program, even if the
    possibility of such damages has been advised against. The entire risk
    as to the quality, the performace, and the fitness of the program for any
    particular purpose lies with the party using the code.

    This code is especially not intended for any clinical or diagnostic use. 
  
 */