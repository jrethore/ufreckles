#include "mex.h"
#include "mat.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *xn,*yn,*zn,*xpix,*ypix,*zpix;
    double *xg,*yg,*zg,*wg;
    
    double dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt, xp, yp, zp;
    double invJ1, invJ2, invJ3, invJ4, invJ5, invJ6, invJ7, invJ8, invJ9;
    double DetJ, dxg, dyg, dzg, res, niter;
    double N_r[8], N_s[8], N_t[8],N[8];
    long dimP, i, j,elt;
    
    mwSize mwDimG[1];
    
    if(nrhs!=7 || nlhs!=4) {
        mexErrMsgTxt("[xg,yg,zg,wg]=GetGaussPointsVoxel(elt,xn,yn,zn,xpix,ypix,zpix)");
    }
    
    elt=mxGetScalar(prhs[0]);
    xn=mxGetData(prhs[1]);
    yn=mxGetData(prhs[2]);
    zn=mxGetData(prhs[3]);
    
    xpix=mxGetData(prhs[4]);
    ypix=mxGetData(prhs[5]);
    zpix=mxGetData(prhs[6]);
    
    dimP=mxGetDimensions(prhs[4])[0];
    
    mwDimG[0] = dimP;
    
    plhs[0] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    
    xg=mxGetData(plhs[0]);
    yg=mxGetData(plhs[1]);
    zg=mxGetData(plhs[2]);
    wg=mxGetData(plhs[3]);
    
    for(i=0;i<dimP;i++) {

        niter = 0;
        res = 1;
        DetJ = 1;
        while((res > 1.e-12) && (niter < 100000) && (DetJ > 0)) {
            
            niter = niter+1;
            
            dxdr = 0;
            dydr = 0;
            dzdr = 0;
            dxds = 0;
            dyds = 0;
            dzds = 0;
            dxdt = 0;
            dydt = 0;
            dzdt = 0;
            xp = 0;
            yp = 0;
            zp = 0;
            
            switch (elt){
                case 8 :
                    N_r[0]=-0.125*(1-(*(yg+i)))*(1-(*(zg+i)));
                    N_r[1]= 0.125*(1-(*(yg+i)))*(1-(*(zg+i)));
                    N_r[2]= 0.125*(1+(*(yg+i)))*(1-(*(zg+i)));
                    N_r[3]=-0.125*(1+(*(yg+i)))*(1-(*(zg+i)));
                    N_r[4]=-0.125*(1-(*(yg+i)))*(1+(*(zg+i)));
                    N_r[5]= 0.125*(1-(*(yg+i)))*(1+(*(zg+i)));
                    N_r[6]= 0.125*(1+(*(yg+i)))*(1+(*(zg+i)));
                    N_r[7]=-0.125*(1+(*(yg+i)))*(1+(*(zg+i)));
                    
                    N_s[0]=-0.125*(1-(*(xg+i)))*(1-(*(zg+i)));
                    N_s[1]=-0.125*(1+(*(xg+i)))*(1-(*(zg+i)));
                    N_s[2]= 0.125*(1+(*(xg+i)))*(1-(*(zg+i)));
                    N_s[3]= 0.125*(1-(*(xg+i)))*(1-(*(zg+i)));
                    N_s[4]=-0.125*(1-(*(xg+i)))*(1+(*(zg+i)));
                    N_s[5]=-0.125*(1+(*(xg+i)))*(1+(*(zg+i)));
                    N_s[6]= 0.125*(1+(*(xg+i)))*(1+(*(zg+i)));
                    N_s[7]= 0.125*(1-(*(xg+i)))*(1+(*(zg+i)));
                    
                    N_t[0]=-0.125*(1-(*(xg+i)))*(1-(*(yg+i)));
                    N_t[1]=-0.125*(1+(*(xg+i)))*(1-(*(yg+i)));
                    N_t[2]=-0.125*(1+(*(xg+i)))*(1+(*(yg+i)));
                    N_t[3]=-0.125*(1-(*(xg+i)))*(1+(*(yg+i)));
                    N_t[4]= 0.125*(1-(*(xg+i)))*(1-(*(yg+i)));
                    N_t[5]= 0.125*(1+(*(xg+i)))*(1-(*(yg+i)));
                    N_t[6]= 0.125*(1+(*(xg+i)))*(1+(*(yg+i)));
                    N_t[7]= 0.125*(1-(*(xg+i)))*(1+(*(yg+i)));
                    
                    N[0]=0.125*(1-(*(xg+i)))*(1-(*(yg+i)))*(1-(*(zg+i)));
                    N[1]=0.125*(1+(*(xg+i)))*(1-(*(yg+i)))*(1-(*(zg+i)));
                    N[2]=0.125*(1+(*(xg+i)))*(1+(*(yg+i)))*(1-(*(zg+i)));
                    N[3]=0.125*(1-(*(xg+i)))*(1+(*(yg+i)))*(1-(*(zg+i)));
                    N[4]=0.125*(1-(*(xg+i)))*(1-(*(yg+i)))*(1+(*(zg+i)));
                    N[5]=0.125*(1+(*(xg+i)))*(1-(*(yg+i)))*(1+(*(zg+i)));
                    N[6]=0.125*(1+(*(xg+i)))*(1+(*(yg+i)))*(1+(*(zg+i)));
                    N[7]=0.125*(1-(*(xg+i)))*(1+(*(yg+i)))*(1+(*(zg+i)));
                    break;
                    
                case 6:
                    N_r[0]=-0.5*(1-(*(zg+i)));
                    N_r[1]= 0.5*(1-(*(zg+i)));
                    N_r[2]= 0;
                    N_r[3]=-0.5*(1+(*(zg+i)));
                    N_r[4]=0.5*(1+(*(zg+i)));
                    N_r[5]=0;
                    N_s[0]=-0.5*(1-(*(zg+i)));
                    N_s[1]= 0;
                    N_s[2]= 0.5*(1-(*(zg+i)));
                    N_s[3]=-0.5*(1+(*(zg+i)));
                    N_s[4]=0;
                    N_s[5]=0.5*(1+(*(zg+i)));
                    N_t[0]=-0.5*(1-(*(xg+i))-(*(yg+i)));
                    N_t[1]=-0.5*(*(xg+i));
                    N_t[2]=-0.5*(*(yg+i));
                    N_t[3]=0.5*(1-(*(xg+i))-(*(yg+i)));
                    N_t[4]=0.5*(*(xg+i));
                    N_t[5]=0.5*(*(yg+i));
                    N[0]=0.5*(1-(*(xg+i))-(*(yg+i)))*(1-(*(zg+i)));
                    N[1]=0.5*(*(xg+i))*(1-(*(zg+i)));
                    N[2]=0.5*(*(yg+i))*(1-(*(zg+i)));
                    N[3]=0.5*(1-(*(xg+i))-(*(yg+i)))*(1+(*(zg+i)));
                    N[4]=0.5*(*(xg+i))*(1+(*(zg+i)));
                    N[5]=0.5*(*(yg+i))*(1+(*(zg+i)));
                    break;
                case 4:
                    N_r[0]=-1;N_r[1]= 1;N_r[2]= 0;N_r[3]=0;
                    N_s[0]=-1;N_s[1]= 0;N_s[2]= 1;N_s[3]=0;
                    N_t[0]=-1;N_t[1]= 0;N_t[2]= 0;N_t[3]=1;
                    N[0]=1-(*(xg+i))-(*(yg+i))-(*(zg+i));
                    N[1]=(*(xg+i));
                    N[2]=(*(yg+i));
                    N[3]=(*(zg+i));
                    break;
                default :
                    break;
                    
            }
            
            
            for(j=0;j<elt;j++) {
                dxdr = dxdr + N_r[j]*(*(xn+j));
                dydr = dydr + N_r[j]*(*(yn+j));
                dzdr = dzdr + N_r[j]*(*(zn+j));
                dxds = dxds + N_s[j]*(*(xn+j));
                dyds = dyds + N_s[j]*(*(yn+j));
                dzds = dzds + N_s[j]*(*(zn+j));
                dxdt = dxdt + N_t[j]*(*(xn+j));
                dydt = dydt + N_t[j]*(*(yn+j));
                dzdt = dzdt + N_t[j]*(*(zn+j));
                xp   = xp   + N[j]*(*(xn+j));
                yp   = yp   + N[j]*(*(yn+j));
                zp   = zp   + N[j]*(*(zn+j));
            }
            
            
            DetJ = dxdr*dyds*dzdt + dxdt*dydr*dzds + dxds*dydt*dzdr - dxdt*dyds*dzdr - dxdr*dydt*dzds - dxds*dydr*dzdt;
            
            if(DetJ > 0) {
                invJ1 =  (dyds*dzdt - dydt*dzds)/DetJ;
                invJ4 = -(dydr*dzdt - dydt*dzdr)/DetJ;
                invJ7 =  (dydr*dzds - dyds*dzdr)/DetJ;
                
                invJ2 = -(dxds*dzdt - dxdt*dzds)/DetJ;
                invJ5 =  (dxdr*dzdt - dxdt*dzdr)/DetJ;
                invJ8 = -(dxdr*dzds - dxds*dzdr)/DetJ;
                
                invJ3 =  (dxds*dydt - dxdt*dyds)/DetJ;
                invJ6 = -(dxdr*dydt - dxdt*dydr)/DetJ;
                invJ9 =  (dxdr*dyds - dxds*dydr)/DetJ;
                
                dxg = invJ1*((*(xpix+i))-xp)+invJ2*((*(ypix+i))-yp)+invJ3*((*(zpix+i))-zp);
                dyg = invJ4*((*(xpix+i))-xp)+invJ5*((*(ypix+i))-yp)+invJ6*((*(zpix+i))-zp);
                dzg = invJ7*((*(xpix+i))-xp)+invJ8*((*(ypix+i))-yp)+invJ9*((*(zpix+i))-zp);
                
                res = dxg*dxg+dyg*dyg+dzg*dzg;
                *(xg+i) = *(xg+i) + dxg;
                *(yg+i) = *(yg+i) + dyg;
                *(zg+i) = *(zg+i) + dzg;
            } else {
                res=0;
                *(xg+i) = 2;
                *(yg+i) = 2;
                *(zg+i) = 2;
                *(wg+i) = 0;
            }
            
        }
        
        switch (elt){
            case 8 :
                if( (fabs(*(xg+i))<1) && (fabs(*(yg+i))<1) && (fabs(*(zg+i))<1) ) {
                    *(wg+i) = 1;
                } else {
                    *(wg+i) = 0;
                    *(xg+i) = 0;
                    *(yg+i) = 0;
                    *(zg+i) = 0;
                }
                break;
                
            case 6:
                break;
                if( ((*(xg+i))>0) && ((*(yg+i))>0) &&  (fabs(*(zg+i))<1) && ((1-(*(xg+i))-(*(yg+i)))>0)  ) {
                    *(wg+i) = 1;
                } else {
                    *(wg+i) = 0;
                    *(xg+i) = 0;
                    *(yg+i) = 0;
                    *(zg+i) = 0;
                }
                
            case 4:
                if( ((*(xg+i))>0) && ((*(yg+i))>0) &&  ((*(zg+i))>0) && ((1-(*(xg+i))-(*(yg+i))-(*(zg+i)))>0)  ) {
                    *(wg+i) = 1;
                } else {
                    *(wg+i) = 0;
                    *(xg+i) = 0;
                    *(yg+i) = 0;
                    *(zg+i) = 0;
                }
                
                break;
            default :
                break;
                
        }

    }
    
}