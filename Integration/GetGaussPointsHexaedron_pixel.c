#include "mex.h"
#include "mat.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double *xn,*yn,*zn,*xpix,*ypix,*zpix;
    double *xg,*yg,*zg,*wg;
    
    double dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt, xp, yp, zp;
    double invJ1, invJ2, invJ3, invJ4, invJ5, invJ6, invJ7, invJ8, invJ9;
    double DetJ, dxg, dyg, dzg, res;
    long N_a[8], N_b[8], N_c[8], dimP, i, j;
    
    mwSize mwDimG[1];
    
    if(nrhs!=6 || nlhs!=4) {
        mexErrMsgTxt("[xg,yg,zg,wg]=func_test_v3(xn,yn,zn,xpix,ypix,zpix)");
    }

    xn=mxGetData(prhs[0]);
    yn=mxGetData(prhs[1]);
    zn=mxGetData(prhs[2]);
    
    xpix=mxGetData(prhs[3]);
    ypix=mxGetData(prhs[4]);
    zpix=mxGetData(prhs[5]);
    
    dimP=mxGetDimensions(prhs[3])[0];
    
    mwDimG[0] = dimP;
    
    plhs[0] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(1, mwDimG, mxDOUBLE_CLASS, mxREAL);

    xg=mxGetData(plhs[0]);
    yg=mxGetData(plhs[1]);
    zg=mxGetData(plhs[2]);
    wg=mxGetData(plhs[3]);
    
    
    N_a[0]=-1;N_a[1]= 1;N_a[2]= 1;N_a[3]=-1;N_a[4]=-1;N_a[5]= 1;N_a[6]= 1;N_a[7]=-1; 
    N_b[0]=-1;N_b[1]=-1;N_b[2]= 1;N_b[3]= 1;N_b[4]=-1;N_b[5]=-1;N_b[6]= 1;N_b[7]= 1; 
    N_c[0]=-1;N_c[1]=-1;N_c[2]=-1;N_c[3]=-1;N_c[4]= 1;N_c[5]= 1;N_c[6]= 1;N_c[7]= 1;
    
    res = 1;
    
    while(res > 1.e-6) {
        
        res = 0;
    
        for(i=0;i<dimP;i++) {
            
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

            for(j=0;j<8;j++) {
                dxdr = dxdr + 0.125*(0+N_a[j]          )*(1+N_b[j]*(*(yg+i)))*(1+N_c[j]*(*(zg+i)))*(*(xn+j));
                dydr = dydr + 0.125*(0+N_a[j]          )*(1+N_b[j]*(*(yg+i)))*(1+N_c[j]*(*(zg+i)))*(*(yn+j));
                dzdr = dzdr + 0.125*(0+N_a[j]          )*(1+N_b[j]*(*(yg+i)))*(1+N_c[j]*(*(zg+i)))*(*(zn+j));
                dxds = dxds + 0.125*(1+N_a[j]*(*(xg+i)))*(0+N_b[j]          )*(1+N_c[j]*(*(zg+i)))*(*(xn+j));
                dyds = dyds + 0.125*(1+N_a[j]*(*(xg+i)))*(0+N_b[j]          )*(1+N_c[j]*(*(zg+i)))*(*(yn+j));
                dzds = dzds + 0.125*(1+N_a[j]*(*(xg+i)))*(0+N_b[j]          )*(1+N_c[j]*(*(zg+i)))*(*(zn+j));
                dxdt = dxdt + 0.125*(1+N_a[j]*(*(xg+i)))*(1+N_b[j]*(*(yg+i)))*(0+N_c[j]          )*(*(xn+j));
                dydt = dydt + 0.125*(1+N_a[j]*(*(xg+i)))*(1+N_b[j]*(*(yg+i)))*(0+N_c[j]          )*(*(yn+j));
                dzdt = dzdt + 0.125*(1+N_a[j]*(*(xg+i)))*(1+N_b[j]*(*(yg+i)))*(0+N_c[j]          )*(*(zn+j));
                xp   = xp   + 0.125*(1+N_a[j]*(*(xg+i)))*(1+N_b[j]*(*(yg+i)))*(1+N_c[j]*(*(zg+i)))*(*(xn+j));
                yp   = yp   + 0.125*(1+N_a[j]*(*(xg+i)))*(1+N_b[j]*(*(yg+i)))*(1+N_c[j]*(*(zg+i)))*(*(yn+j));
                zp   = zp   + 0.125*(1+N_a[j]*(*(xg+i)))*(1+N_b[j]*(*(yg+i)))*(1+N_c[j]*(*(zg+i)))*(*(zn+j));
            }
        
            DetJ = dxdr*dyds*dzdt + dxdt*dydr*dzds + dxds*dydt*dzdr - dxdt*dyds*dzdr - dxdr*dydt*dzds - dxds*dydr*dzdt;

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

            res = res + dxg*dxg+dyg*dyg+dzg*dzg;

            *(xg+i) = *(xg+i) + dxg;
            *(yg+i) = *(yg+i) + dyg;
            *(zg+i) = *(zg+i) + dzg;
        
        }
    
    }
    
    for(i=0;i<dimP;i++) {
       if( (abs(*(xg+i))<1) && (abs(*(yg+i))<1) && (abs(*(zg+i))<1) ) {
            *(wg+i) = 1;  
        }
    }
    
}   