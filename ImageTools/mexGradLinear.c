#include "mex.h"
#include "mat.h"
#include "math.h"

void mexGradLinear(double *Upp,double *Vpp,double *im1x,double *im1y, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc)
{
    int ii,jj,IX1,IY1,ind,count=0;
    double Xi,Yi;
    for (jj=0;jj<n;jj++){
        for (ii=0;ii<m;ii++){
            Xi=*(Upp+count);
            Yi=*(Vpp+count);
            IX1 = floor(Xi);
            IY1 = floor(Yi);
            ind=IX1+(IY1-1)*mc-1;
            if ((IX1-1<0)||(IX1-1>mc-1)||(IY1-1<0)||(IY1-1>nc-1)){
                *(im1x+count)   = 0.;
                *(im1y+count)   = 0.;
            }
            else {
                if (IX1-1==0){
                    *(im1x+count)   = (*(im0+ind+1)-*(im0+ind));
                }
                else if (IX1-1==mc-1){
                    *(im1x+count)   = (*(im0+ind)-*(im0+ind-1));
                }
                else {
                    *(im1x+count)   = 0.5*(*(im0+ind+1)-*(im0+ind-1));
                }
                if (IY1-1==0){
                    *(im1y+count)   = (*(im0+ind+mc)-*(im0+ind));
                }
                else if (IY1-1==nc-1){
                    *(im1y+count)   = (*(im0+ind)-*(im0+ind-mc));
                }
                else {
                    *(im1y+count)   = 0.5*(*(im0+ind+mc)-*(im0+ind-mc));
                }
            }
            count++;
        }
        
    }
    
    
}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{

    double *im0,*im1x,*im1y;
    double *up,*vp;
    mwSize mrows0,ncols0;
    mwSize mrows1,ncols1;
    mwSize dims[2];


   if(nrhs<3)
        mexErrMsgTxt("At least three inputs required.");

    if(nlhs!=2)
        mexErrMsgTxt("Only two outputs required.");

        up= mxGetData(prhs[0]);
    vp= mxGetData(prhs[1]);
im0= mxGetData(prhs[2]);
    mrows1 = mxGetM(prhs[0]);
    ncols1 = mxGetN(prhs[0]);
    mrows0 = mxGetM(prhs[2]);
    ncols0 = mxGetN(prhs[2]);
     dims[0]=mrows1;dims[1]=ncols1;;

    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);    
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);    
im1x= mxGetData(plhs[0]);
im1y= mxGetData(plhs[1]);

    mexGradLinear(up,vp,im1x,im1y,mrows1,ncols1,im0,mrows0,ncols0);


}
