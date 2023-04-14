#include "mex.h"
#include "mat.h"

void FDGradient(double *im,double *gx, mwSize m, mwSize n)
{
    int i,j,count=0;
    for (j=0;j<n;j++){
        
        *(gx+count)=*(im+count+1)-*(im+count);
        
        count++;
        for (i=1;i<m-1;i++){
            *(gx+count)=0.5*(*(im+count+1)-*(im+count-1));
            count++;
            
        }
        *(gx+count)=*(im+count)-*(im+count-1);
        count++;
    }
    
}

void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *im,*gx;
    mwSize mrows,ncols;
     mwSize  nsubs;
     if(nrhs!=1)
        mexErrMsgTxt("One inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");
    
    im= mxGetPr(prhs[0]);
     nsubs=mxGetNumberOfDimensions(prhs[0]);
     mrows = mxGetM(prhs[0]);
         ncols = mxGetN(prhs[0]);
    ncols = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);    
gx=mxGetPr(plhs[0]);
    FDGradient(im,gx,mrows,ncols);

}
