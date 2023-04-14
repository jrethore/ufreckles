#include "mex.h"
#include "mat.h"
#include "math.h"
void GaussPointsPixels(double *xg, double *yg, double *wg, int elt, double *xpix,double *ypix, double *xn, double *yn, mwSize dimP)
{
    double dxdr, dydr, dxds, dyds, xp, yp;
    double invJ1, invJ2, invJ3, invJ4;
    double DetJ, dxg, dyg, res, niter;
    double N_r[4], N_s[4], N_t[4],N[4];
    long i, j;
    for(i=0;i<dimP;i++) {
        
        niter = 0;
        res = 1;
        DetJ = 1;
        
        while((res > 1.e-12) && (niter < 1000) && (DetJ > 0)) {
            
            niter = niter+1;
            
            dxdr = 0;
            dydr = 0;
            dxds = 0;
            dyds = 0;
            xp = 0;
            yp = 0;
            
            switch (elt){
                case 4:
                    N_r[0]=-0.25*(1-(*(yg+i)));
                    N_r[1]= 0.25*(1-(*(yg+i)));
                    N_r[2]= 0.25*(1+(*(yg+i)));
                    N_r[3]=-0.25*(1+(*(yg+i)));
                    
                    N_s[0]=-0.25*(1-(*(xg+i)));
                    N_s[1]=-0.25*(1+(*(xg+i)));
                    N_s[2]= 0.25*(1+(*(xg+i)));
                    N_s[3]= 0.25*(1-(*(xg+i)));
                    
                    
                    N[0]=0.25*(1-(*(xg+i)))*(1-(*(yg+i)));
                    N[1]=0.25*(1+(*(xg+i)))*(1-(*(yg+i)));
                    N[2]=0.25*(1+(*(xg+i)))*(1+(*(yg+i)));
                    N[3]=0.25*(1-(*(xg+i)))*(1+(*(yg+i)));
                    break;
                    
                case 3 :
                    N_r[0]=-1;N_r[1]= 1;N_r[2]= 0;
                    N_s[0]=-1;N_s[1]= 0;N_s[2]= 1;
                    N[0]=1-(*(xg+i))-(*(yg+i));
                    N[1]=(*(xg+i));
                    N[2]=(*(yg+i));
                    break;
                default :
                    break;
            }
            
            
            for(j=0;j<elt;j++) {
                dxdr = dxdr + N_r[j]*(*(xn+j));
                dydr = dydr + N_r[j]*(*(yn+j));
                dxds = dxds + N_s[j]*(*(xn+j));
                dyds = dyds + N_s[j]*(*(yn+j));
                xp   = xp   + N[j]*(*(xn+j));
                yp   = yp   + N[j]*(*(yn+j));
            }
            
            
            DetJ = dxdr*dyds - dxds*dydr;
            
            
            if(DetJ > 0) {
                invJ1 =  dyds/DetJ;
                invJ2 = -dxds/DetJ;
                invJ3 = -dydr/DetJ;
                invJ4 =  dxdr/DetJ;
                
                dxg = invJ1*((*(xpix+i))-xp)+invJ2*((*(ypix+i))-yp);
                dyg = invJ3*((*(xpix+i))-xp)+invJ4*((*(ypix+i))-yp);
                
                res = dxg*dxg+dyg*dyg;
                
                *(xg+i) = *(xg+i) + dxg;
                *(yg+i) = *(yg+i) + dyg;
            } else {
                res=0;
                *(xg+i) = 2;
                *(yg+i) = 2;
                *(wg+i) = 0;
            }
            
            /*mexPrintf("res = %e\n", res);*/
            
        }
        
        switch (elt){
            case 4:
                if( (fabs(*(xg+i))<1) && (fabs(*(yg+i))<1) ) {
                    *(wg+i) = 1;
                } else {
                    *(wg+i) = 0;
                    *(xg+i) = 0;
                    *(yg+i) = 0;
                }
                break;
                
            case 3 :
                if( ((*(xg+i))>0) && ((*(yg+i))>0) && ((1-(*(xg+i))-(*(yg+i)))>0) ) {
                    *(wg+i) = 1;
                } else {
                    *(wg+i) = 0;
                    *(xg+i) = 0;
                    *(yg+i) = 0;
                }
                break;
            default :
                break;
        }
        
    }
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *xn,*yn,*xpix,*ypix;
    double *xg,*yg,*wg;
    long elt;
    
    mwSize mwDimG[1],dimP;
    
    if(nrhs!=5 || nlhs!=3) {
        mexErrMsgTxt("[xg,yg,wg]=GetGaussPointsPixel(elt,xn,yn,xpix,ypix)");
    }
    
    elt=mxGetScalar(prhs[0]);
    xn=mxGetData(prhs[1]);
    yn=mxGetData(prhs[2]);
    
    xpix=mxGetData(prhs[3]);
    ypix=mxGetData(prhs[4]);
    
    dimP=mxGetDimensions(prhs[3])[0];
    mwDimG[0] = dimP;
    
    
    plhs[0] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    
    xg=mxGetData(plhs[0]);
    yg=mxGetData(plhs[1]);
    wg=mxGetData(plhs[2]);
    
    GaussPointsPixels(xg,yg,wg,elt,xpix,ypix,xn,yn,dimP);
    
}