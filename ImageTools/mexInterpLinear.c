#include "mex.h"
#include "mat.h"
#include "math.h"

void mexInterpLinear(double *Upp,double *Vpp,double *im1, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc)
{
    int ii,jj,IX1,IY1,mcnc,ind,count;
    double Xi,Yi,Xx,Yy,NivGris;
double CR[16];
    mcnc=mc*nc;
    for (jj=0;jj<n;jj++){
        for (ii=0;ii<m;ii++){
            count=ii+(jj)*m;
/*             Xi=*(Upp+count)+ii+1;
             Yi=*(Vpp+count)+jj+1;*/
            Xi=*(Upp+count);
            Yi=*(Vpp+count);
            IX1 = floor(Xi);
            IY1 = floor(Yi);
            Xx=2*(Xi-IX1-0.5);
            Yy=2*(Yi-IY1-0.5);
            ind=IX1+(IY1-1)*mc-1;
                 if ((IX1-1>0)&&(IY1-1>0)&&(IX1-1<mc-1)&&(IY1-1<nc-1)){
       
            NivGris=(*(im0+ind))*0.25*(1-Xx)*(1-Yy)+
            (*(im0+ind+1))*0.25*(1+Xx)*(1-Yy)+
            (*(im0+ind+mc))*0.25*(1-Xx)*(1+Yy)+
            (*(im0+ind+mc+1))*0.25*(1+Xx)*(1+Yy);

            
            *(im1+count)   = NivGris;
                 }
            else {
                *(im1+count)   = -(2^16-1);
            }
        }
        
    }
    
    
}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{

    double *im0,*im1;
    double *up,*vp;
    mwSize mrows0,ncols0;
    mwSize mrows1,ncols1;
    mwSize dims[2];


   if(nrhs<3)
        mexErrMsgTxt("At least three inputs required.");

    if(nlhs!=1)
        mexErrMsgTxt("Only one output required.");

        up= mxGetData(prhs[0]);
    vp= mxGetData(prhs[1]);
im0= mxGetData(prhs[2]);
    mrows1 = mxGetM(prhs[0]);
    ncols1 = mxGetN(prhs[0]);
    mrows0 = mxGetM(prhs[2]);
    ncols0 = mxGetN(prhs[2]);
     dims[0]=mrows1;dims[1]=ncols1;;

    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);    
im1= mxGetData(plhs[0]);

    mexInterpLinear(up,vp,im1,mrows1,ncols1,im0,mrows0,ncols0);


}
