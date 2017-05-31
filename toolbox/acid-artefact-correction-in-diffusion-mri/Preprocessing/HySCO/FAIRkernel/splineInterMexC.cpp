/*
 * (c) Jan Modersitzki and Fabian Gigengack 2011/02/02, see FAIR.2 and FAIRcopyright.m.
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for spline interpolation. See splineInter for details.
 */

#include <mex.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

void splineInter1D(double *Tc,double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative);
void splineInter2D(double *Tc,double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative);
void splineInter3D(double *Tc,double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative);
inline double  b0(int j, double xi);
inline double db0(int j, double xi);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int dim, i;
	if (nrhs<5)
		mexErrMsgTxt("Number of arguments must be 5!");
	
    // Get input
	const double* T     = static_cast<double*>(mxGetData(prhs[0]));
	const double* omega = static_cast<double*>(mxGetData(prhs[1]));
	const double* m     = static_cast<double*>(mxGetData(prhs[2]));
    const double* X     = static_cast<double*>(mxGetData(prhs[3]));
	bool doDerivative   = (bool) *mxGetLogicals(prhs[4]);
    
	// get dimension, h and number of elements N
	dim = mxGetN(prhs[1]) / 2;
	double h[3]; //h[dim];
	for (i=0; i<dim; i++) {
		h[i] = (omega[2*i+1]-omega[2*i]) / m[i];
    }
	const int N = mxGetM(prhs[3])/dim;
    
	//Allocate Tc
	mwSize dims[2]; dims[0] = N; dims[1] = 1;
	plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
	double* Tc = static_cast<double*>(mxGetData(plhs[0]));
	
	double* dT = 0;
    //Allocate dT
	if (doDerivative) {
        dims[1] = dim; 
        plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
		dT = static_cast<double*>(mxGetData(plhs[1]));
    } else {
        dims[0] = 1; dims[1] = 1;
        plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
        dT = static_cast<double*>(mxGetData(plhs[1]));
    }
    
    switch (dim) {
		case 1:
			splineInter1D(Tc,dT,T,omega,m,N,X,h,doDerivative);
			break;
		case 2:
			splineInter2D(Tc,dT,T,omega,m,N,X,h,doDerivative);
			break;
		case 3:
			splineInter3D(Tc,dT,T,omega,m,N,X,h,doDerivative);
			break;
		default:
			break;
	}
}

void splineInter1D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative) {
    int xf, i, j, m0 = m[0];
    double x, p;
    
    #pragma omp parallel for default(shared) private(x, xf, p, i, j)
    for (i=0;i<N;i++) {
        x = (X[i] - omega[0]) / h[0] + .5 - 1; //subtract 1 for indexing purposes in C
        
        //check if it is a valid x
        if (x<=-2 || x>=m0+1)
            continue;
        
        xf = floor(x);
        x  = x - xf;
        
        Tc[i] = 0;
        if (doDerivative){
            dT[i] = 0;
        }
        for (j=-1;j<3;j++) {
            p = (xf+j<0 || xf+j>m0-1)? 0: T[xf+j];
            Tc[i] += p*b0(3-j,x-j);
            if (doDerivative) {
                dT[i] += p*db0(3-j,x-j);
            }
        }
        
        if (doDerivative) {
            dT[i] = dT[i]/h[0];
        }
    }
}

void splineInter2D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative) {
	// Increments in X and Y direction
    int xf, yf, i, j, k, m0 = m[0], m1 = m[1];
	double x, y, p;
    
    //for each value X compute the value in  the spline
    #pragma omp parallel for default(shared) private(x, y, xf, yf, p, i, j, k)
    for (i=0;i<N;i++) {
        x = (X[i]   - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N] - omega[2]) / h[1] + .5 - 1;

        //check if it is a valid x
        if (x<=-2 || y<=-2 || x>=m0+1 || y>=m1+1)
            continue;
        
        xf = floor(x);
        yf = floor(y);
        x  = x - xf;
        y  = y - yf;
        
        Tc[i]   = 0;
        if (doDerivative){
            dT[i]   = 0;
            dT[i+N] = 0;
        }
        for (j=-1;j<3;j++) {
            for (k=-1;k<3;k++) {
                p = (xf+j<0 || xf+j>m0-1 ||
                     yf+k<0 || yf+k>m1-1)? 0: T[xf+j+m0*(yf+k)];
                Tc[i] += p*b0(3-j,x-j)*b0(3-k,y-k);
                if (doDerivative) {
                    dT[i]   += p*db0(3-j,x-j)* b0(3-k,y-k);
                    dT[i+N] += p* b0(3-j,x-j)*db0(3-k,y-k);
                }
            }
        }
        if (doDerivative) {
            dT[i]   = dT[i]/h[0];
            dT[i+N] = dT[i+N]/h[1];
        }
    }	
}

void splineInter3D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative) {
    int xf, yf, zf, i, j, k, l;
    int i1 = 1, m0 = m[0], m1 = m[1], m2 = m[2], m01 = m0*m1;
    double x, y, z, p;
    
	// Increments in X,Y and Z direction
    #pragma omp parallel for default(shared) private(x, y, z, xf, yf, zf, p, i, j, k, l)
    for (i=0;i<N;i++) {
        x = (X[i]     - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N]   - omega[2]) / h[1] + .5 - 1;
        z = (X[i+2*N] - omega[4]) / h[2] + .5 - 1;
		
		xf = floor(x);
        yf = floor(y);
        zf = floor(z);
        //check if it is a valid x
        if (x<=-2 || y<=-2 || z<=-2 || x>=m[0]+1 || y>=m1+1|| z>=m2+1)
            continue;
        
        // Extract remainder
        x = x - xf;
        y = y - yf;
        z = z - zf;
        
        Tc[i]     = 0;
        if(doDerivative)
        {
            dT[i]     = 0;
            dT[i+N]   = 0;
            dT[i+2*N] = 0;
        }
		for (j=-1;j<3;j++) {
            for (k=-1;k<3;k++) {
				for (l=-1;l<3;l++) {
                    p = (xf+j<0 || xf+j>m0-1 ||
                         yf+k<0 || yf+k>m1-1 ||
                         zf+l<0 || zf+l>m2-1)? 0: T[xf+j+m0*(yf+k)+m01*(zf+l)];
					Tc[i] += p*b0(3-j,x-j)*b0(3-k,y-k)*b0(3-l,z-l);
					if (doDerivative) {
                        dT[i]     += p*db0(3-j,x-j)* b0(3-k,y-k)* b0(3-l,z-l);
						dT[i+N]   += p* b0(3-j,x-j)*db0(3-k,y-k)* b0(3-l,z-l);
						dT[i+2*N] += p* b0(3-j,x-j)* b0(3-k,y-k)*db0(3-l,z-l);
					}
				}
            }
        }
        if (doDerivative) {
            dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
        }
    }
}


inline double b0(int j, double xi) {
    switch (j) {
        case 1:
            return (2+xi) * (2+xi) * (2+xi); //  (2+xi).^3;
        case 2:
            return -(3*xi+6)*xi*xi+4;        // -(3*xi+6).*xi.^2+4;
        case 3:
            return  (3*xi-6)*xi*xi+4;        //  (3*xi-6).*xi.^2+4;
        case 4:
            return (2-xi) * (2-xi) * (2-xi); //  (2-xi).^3;
    }
    return 0;
}

inline double db0(int j, double xi) {
    switch (j) {
        case 1:
            return  3*(2+xi)*(2+xi);         //  3*(2+xi).^2;
        case 2:
            return -(9*xi+12)*xi;            // -(9*xi+12).*xi;
        case 3:
            return  (9*xi-12)*xi;            //  (9*xi-12).*xi;
        case 4:  
            return -3*(2-xi)*(2-xi);         // -3*(2-xi).^2;
    }
    return 0;
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