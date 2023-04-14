#include "mex.h"
#include "mat.h"
#include "math.h"

void mexInterpSpline(double *Upp,double *Vpp,double *Wpp,double *im1, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc, mwSize oc)
{
    int ii,jj,IX1,IY1,IZ1,mcnc,ind,count=0;
    double Xi,Yi,Zi,Xx,Yy,Zz;
    double NivGris;
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
            Xx=2*(Xi-IX1-0.5);
            Yy=2*(Yi-IY1-0.5);
            Zz=2*(Zi-IZ1-0.5);
            if ((IX1-1>0)&&(IY1-1>0)&&(IZ1-1>0)&&(IX1-1<mc-1)&&(IY1-1<nc-1)&&(IZ1-1<oc-1)){
                
                ind=IX1-1+((IY1-1)+(IZ1-1)*nc)*mc;
                NivGris=(*(im0+ind+ijk))*0.125*(1-Xx)*(1-Yy)*(1-Zz)+
                        (*(im0+ind+ijk+1))*0.125*(1+Xx)*(1-Yy)*(1-Zz)+
                        (*(im0+ind+ijpk))*0.125*(1-Xx)*(1+Yy)*(1-Zz)+
                        (*(im0+ind+ijpk+1))*0.125*(1+Xx)*(1+Yy)*(1-Zz)+
                        (*(im0+ind+ijkp))*0.125*(1-Xx)*(1-Yy)*(1+Zz)+
                        (*(im0+ind+ijkp+1))*0.125*(1+Xx)*(1-Yy)*(1+Zz)+
                        (*(im0+ind+ijpkp))*0.125*(1-Xx)*(1+Yy)*(1+Zz)+
                        (*(im0+ind+ijpkp+1))*0.125*(1+Xx)*(1+Yy)*(1+Zz);
                
                *(im1+count)   = NivGris;
            }
            else {
                *(im1+count)   = 0;
            }
            count++;
        }
        
    }
    
}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{

    double *im0,*im1;
    double *up,*vp,*wp;
    mwSize mrows0,ncols0,pstacks0;
    mwSize mrows1,ncols1;
    mwSize dims[2];

   if(nrhs<4)
        mexErrMsgTxt("At least four inputs required.");

    if(nlhs!=1)
        mexErrMsgTxt("Only one output required.");

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
im1= mxGetData(plhs[0]);

    mexInterpSpline(up,vp,wp,im1,mrows1,ncols1,im0,mrows0,ncols0,pstacks0);


}
