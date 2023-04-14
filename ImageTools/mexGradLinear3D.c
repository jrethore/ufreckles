#include "mex.h"
#include "mat.h"
#include "math.h"

void mexGradSpline(double *Upp,double *Vpp,double *Wpp,double *gx,double *gy,double *gz, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc, mwSize oc)
{
    
    int ii,jj,IX1,IY1,IZ1,mcnc,ind,count=0;
    double Xi,Yi,Zi,Xx,Yy,Zz;
    double dx,dy,dz;
    int ijk,ijkp;
    int ijpk,ijpkp;
    mcnc=mc*nc;
    ijk = 0;
    ijkp = mcnc;
    ijpk = mc;
    ijpkp = mc+mcnc;
    
    
    for (jj=0;jj<n;jj++){
        for (ii=0;ii<m;ii++){
            Xi=*(Upp+count);
            Yi=*(Vpp+count);
            Zi=*(Wpp+count);
            IX1 = floor(Xi);
            IY1 = floor(Yi);
            IZ1 = floor(Zi);
            ind=IX1-1+((IY1-1)+(IZ1-1)*nc)*mc;
            
            if ((IX1-1<0)||(IX1-1>mc-1)||(IY1-1<0)||(IY1-1>nc-1)||(IZ1-1<0)||(IZ1-1>oc-1)){
                dx=0.;
                dy=0.;
                dz=0.;
            }
            else {
                if (IX1-1==0){
                    dx   = (*(im0+ind+1)-*(im0+ind));
                }
                else if (IX1-1==mc-1){
                    dx = (*(im0+ind)-*(im0+ind-1));
                }
                else {
                    dx  = 0.5*(*(im0+ind+1)-*(im0+ind-1));
                }
                
                if (IY1-1==0){
                    dy   = (*(im0+ind+mc)-*(im0+ind));
                }
                else if (IY1-1==nc-1){
                    dy = (*(im0+ind)-*(im0+ind-mc));
                }
                else {
                    dy  = 0.5*(*(im0+ind+mc)-*(im0+ind-mc));
                }
                if (IZ1-1==0){
                    dz   = (*(im0+ind+mcnc)-*(im0+ind));
                }
                else if (IZ1-1==oc-1){
                    dz = (*(im0+ind)-*(im0+ind-mcnc));
                }
                else {
                    dz  = 0.5*(*(im0+ind+mcnc)-*(im0+ind-mcnc));
                }
                
            }
            *(gx+count)   = dx;
            *(gy+count)   = dy;
            *(gz+count)   = dz;
            
            count++;
        }
        
    }
    
}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    
    double *im0,*gx,*gy,*gz;
    double *up,*vp,*wp;
    mwSize mrows0,ncols0,pstacks0;
    mwSize mrows1,ncols1;
    mwSize dims[2];
    
    if(nrhs<4)
        mexErrMsgTxt("At least four inputs required.");
    
    if(nlhs!=3)
        mexErrMsgTxt("Three outputs required.");
    
    up= mxGetData(prhs[0]);
    vp= mxGetData(prhs[1]);
    wp= mxGetData(prhs[2]);
    im0= mxGetData(prhs[3]);
    mrows1 = mxGetM(prhs[0]);
    ncols1 = mxGetN(prhs[0]);
    mrows0 = mxGetM(prhs[3]);
    ncols0 = ((mxGetDimensions(prhs[3]))[1]);
    pstacks0 = ((mxGetDimensions(prhs[3]))[2]);
    dims[0]=mrows1;
    dims[1]=ncols1;
    
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    gx= mxGetData(plhs[0]);
    gy= mxGetData(plhs[1]);
    gz= mxGetData(plhs[2]);
    
    mexGradSpline(up,vp,wp,gx,gy,gz,mrows1,ncols1,im0,mrows0,ncols0,pstacks0);
    
    
}
