#include "mex.h"
#include "mat.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *xg,*yg,*FF;
    long dimP, i, j,elt;
    int count=0;
    mwSize mwDimG[2];
    
    if(nrhs!=3 || nlhs!=1) {
        mexErrMsgTxt("[FF]=GetFiniteElementShapeFunctions(elt,xg,yg)");
    }
    
    elt=mxGetScalar(prhs[0]);
    xg=mxGetData(prhs[1]);
    yg=mxGetData(prhs[2]);
    
    
    dimP=mxGetDimensions(prhs[1])[0];
    
    mwDimG[0] = dimP;
    mwDimG[1] = 4;
    plhs[0] = mxCreateNumericArray(2, mwDimG, mxDOUBLE_CLASS, mxREAL);
    FF= mxGetData(plhs[0]);

    
    
    
    switch (elt){
        case 4 :
            
                for(i=0;i<dimP;i++) {
                    *(FF + count)=0.25*(1-(*(xg+i)))*(1-(*(yg+i)));
                    count++;
                }
                for(i=0;i<dimP;i++) {
                    *(FF + count)=0.25*(1+(*(xg+i)))*(1-(*(yg+i)));
                    count++;
                }
                for(i=0;i<dimP;i++) {
                    *(FF + count)=0.25*(1+(*(xg+i)))*(1+(*(yg+i)));
                    count++;
                }
                for(i=0;i<dimP;i++) {
                    *(FF + count)=0.25*(1-(*(xg+i)))*(1+(*(yg+i)));
                    count++;
                }
                
            
            break;
            
        case 3:
                for(i=0;i<dimP;i++) {
                    *(FF + count)=1-(*(xg+i))-(*(yg+i));
                    count++;
                }
                for(i=0;i<dimP;i++) {
                    *(FF + count)=(*(xg+i));
                    count++;
                }
                for(i=0;i<dimP;i++) {
                    *(FF + count)=(*(yg+i));
                    count++;
                }
                
            break;
        default :
            break;
            
    }
    
    
}