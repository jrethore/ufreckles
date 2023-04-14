#include "mex.h"
#include "mat.h"
#include "math.h"

void mexPhidf2D(double *phidf,double *N,double *im, mwSize m, mwSize n,mwSize m1,mwSize n1,int step)
{
    
    int i,j,l,count=0;
    int pix;
    double gx,gy;
    
    for (j=0;j<n;j+=step){
        for (i=0;i<m;i+=step){
            pix=i+j*m;
            if (i==0) gx=*(im+pix+1)-*(im+pix);
            else if (i==m-1) gx=*(im+pix)-*(im+pix-1);
            else gx=0.5*(*(im+pix+1)-*(im+pix-1));
            if (j==0) gy=*(im+pix+m)-*(im+pix);
            else if (j==n-1) gy=*(im+pix)-*(im+pix-m);
            else gy=0.5*(*(im+pix+m)-*(im+pix-m));
            
            
            for (l=0;l<n1;l++){
                *(phidf+count+l*m1)=*(N+count+l*m1)*gx;
                *(phidf+count+l*m1+n1*m1)=*(N+count+l*m1)*gy;
            }
            count++;
        }
    }
    
}






void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    double *im0;
    double *N,*phidf;
    mwSize mrows0,ncols0;
    mwSize mrows1,ncols1;
    mwSize dims[2];
    int step;
    if (nrhs==2){
        step=1;
    } else {
        step=mxGetScalar(prhs[2]);
    }
    
    if(nrhs<2)
        mexErrMsgTxt("At least two inputs required.");
    
    if(nlhs!=1)
        mexErrMsgTxt("Only one output required.");
    
    im0= mxGetData(prhs[0]);
    N= mxGetData(prhs[1]);
    mrows0 = mxGetM(prhs[0]);
    ncols0 = ((mxGetDimensions(prhs[0]))[1]);
    mrows1 = mxGetM(prhs[1]);
    ncols1 = mxGetN(prhs[1]);
    dims[0]=mrows1;
    dims[1]=2*ncols1;
    
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    phidf= mxGetData(plhs[0]);
    
    mexPhidf2D(phidf,N,im0,mrows0,ncols0,mrows1,ncols1,step);
    
    
}
