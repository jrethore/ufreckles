#include "mex.h"
#include "mat.h"
#include "math.h"

void mexInterpSpline(double *Upp,double *Vpp,double *Wpp,double *im1, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc, mwSize oc)
{
    int ii,jj,IX1,IY1,IZ1,mcnc,ind,count=0;
    double Xi,Yi,Zi,Xx,Yy,Zz;
    double P0X,P1X,P2X,P3X;
    double P0Y,P1Y,P2Y,P3Y;
    double P0Z,P1Z,P2Z,P3Z;
    double NivGris0,NivGris1,NivGris2,NivGris3;
    int ijmkm,ijmk,ijmkp,ijmkpp;
    int ijkm,ijk,ijkp,ijkpp;
    int ijpkm,ijpk,ijpkp,ijpkpp;
    int ijppkm,ijppk,ijppkp,ijppkpp;
    mcnc=mc*nc;

    ijmkm = -mc-mcnc;
    ijmk = -mc;
    ijmkp = -mc+mcnc;
    ijmkpp = -mc+2*mcnc;
    ijkm = -mcnc;
    ijk = 0;
    ijkp = mcnc;
    ijkpp = 2*mcnc;
    ijpkm = mc-mcnc;
    ijpk = mc;
    ijpkp = mc+mcnc;
    ijpkpp = mc+2*mcnc;
    ijppkm = 2*mc-mcnc;
    ijppk = 2*mc;
    ijppkp = 2*mc+mcnc;
    ijppkpp = 2*mc+2*mcnc;
    
    
    for (jj=0;jj<n;jj++){
        for (ii=0;ii<m;ii++){
/*             Xi=*(Upp+count)+ii+1;
             Yi=*(Vpp+count)+jj+1;*/
            Xi=*(Upp+count);
            Yi=*(Vpp+count);
            Zi=*(Wpp+count);
            IX1 = floor(Xi);
            IY1 = floor(Yi);
            IZ1 = floor(Zi);
            Xx=Xi-IX1;
            Yy=Yi-IY1;
            Zz=Zi-IZ1;
            ind=IX1+((IY1-1)+(IZ1-1)*nc)*mc-1;
            
            P0X=0.2*Xx*(-7/3+Xx*(4-5*Xx/3));
            P1X=1+0.2*Xx*(-1+Xx*(-9+5*Xx));
            P2X=0.2*Xx*(4+Xx*(6-5*Xx));
            P3X=0.2*Xx*(-2/3+Xx*(-1+5/3*Xx));
            
            P0Y=0.2*Yy*(-7/3+Yy*(4-5*Yy/3));
            P1Y=1+0.2*Yy*(-1+Yy*(-9+5*Yy));
            P2Y=0.2*Yy*(4+Yy*(6-5*Yy));
            P3Y=0.2*Yy*(-2/3+Yy*(-1+5/3*Yy));

            P0Z=0.2*Zz*(-7/3+Zz*(4-5*Zz/3));
            P1Z=1+0.2*Zz*(-1+Zz*(-9+5*Zz));
            P2Z=0.2*Zz*(4+Zz*(6-5*Zz));
            P3Z=0.2*Zz*(-2/3+Zz*(-1+5/3*Zz));
            
                      
            
            NivGris0=P0X*(
                           P0Y*(
                                P0Z*(*(im0+ind-1+ijmkm))+P1Z*(*(im0+ind-1+ijmk))+P2Z*(*(im0+ind-1+ijmkp))+P3Z*(*(im0+ind-1+ijmkpp))
                               )
                           +P1Y*(
                                P0Z*(*(im0+ind-1+ijkm))+P1Z*(*(im0+ind-1+ijk))+P2Z*(*(im0+ind-1+ijkp))+P3Z*(*(im0+ind-1+ijkpp))
                               )
                           +P2Y*(
                                P0Z*(*(im0+ind-1+ijpkm))+P1Z*(*(im0+ind-1+ijpk))+P2Z*(*(im0+ind-1+ijpkp))+P3Z*(*(im0+ind-1+ijpkpp))
                               )
                           +P3Y*(
                                P0Z*(*(im0+ind-1+ijppkm))+P1Z*(*(im0+ind-1+ijppk))+P2Z*(*(im0+ind-1+ijppkp))+P3Z*(*(im0+ind-1+ijppkpp))
                               )                           
                          );
                      
            NivGris1=P1X*(
                           P0Y*(
                                P0Z*(*(im0+ind+ijmkm))+P1Z*(*(im0+ind+ijmk))+P2Z*(*(im0+ind+ijmkp))+P3Z*(*(im0+ind+ijmkpp))
                               )
                           +P1Y*(
                                P0Z*(*(im0+ind+ijkm))+P1Z*(*(im0+ind+ijk))+P2Z*(*(im0+ind+ijkp))+P3Z*(*(im0+ind+ijkpp))
                               )
                           +P2Y*(
                                P0Z*(*(im0+ind+ijpkm))+P1Z*(*(im0+ind+ijpk))+P2Z*(*(im0+ind+ijpkp))+P3Z*(*(im0+ind+ijpkpp))
                               )
                           +P3Y*(
                                P0Z*(*(im0+ind+ijppkm))+P1Z*(*(im0+ind+ijppk))+P2Z*(*(im0+ind+ijppkp))+P3Z*(*(im0+ind+ijppkpp))
                               )                           
                          );
                                
            NivGris2=P2X*(
                           P0Y*(
                                P0Z*(*(im0+ind+1+ijmkm))+P1Z*(*(im0+ind+1+ijmk))+P2Z*(*(im0+ind+1+ijmkp))+P3Z*(*(im0+ind+1+ijmkpp))
                               )
                           +P1Y*(
                                P0Z*(*(im0+ind+1+ijkm))+P1Z*(*(im0+ind+1+ijk))+P2Z*(*(im0+ind+1+ijkp))+P3Z*(*(im0+ind+1+ijkpp))
                               )
                           +P2Y*(
                                P0Z*(*(im0+ind+1+ijpkm))+P1Z*(*(im0+ind+1+ijpk))+P2Z*(*(im0+ind+1+ijpkp))+P3Z*(*(im0+ind+1+ijpkpp))
                               )
                           +P3Y*(
                                P0Z*(*(im0+ind+1+ijppkm))+P1Z*(*(im0+ind+1+ijppk))+P2Z*(*(im0+ind+1+ijppkp))+P3Z*(*(im0+ind+1+ijppkpp))
                               )                           
                          );
                                
            NivGris3=P3X*(
                           P0Y*(
                                P0Z*(*(im0+ind+2+ijmkm))+P1Z*(*(im0+ind+2+ijmk))+P2Z*(*(im0+ind+2+ijmkp))+P3Z*(*(im0+ind+2+ijmkpp))
                               )
                           +P1Y*(
                                P0Z*(*(im0+ind+2+ijkm))+P1Z*(*(im0+ind+2+ijk))+P2Z*(*(im0+ind+2+ijkp))+P3Z*(*(im0+ind+2+ijkpp))
                               )
                           +P2Y*(
                                P0Z*(*(im0+ind+2+ijpkm))+P1Z*(*(im0+ind+2+ijpk))+P2Z*(*(im0+ind+2+ijpkp))+P3Z*(*(im0+ind+2+ijpkpp))
                               )
                           +P3Y*(
                                P0Z*(*(im0+ind+2+ijppkm))+P1Z*(*(im0+ind+2+ijppk))+P2Z*(*(im0+ind+2+ijppkp))+P3Z*(*(im0+ind+2+ijppkpp))
                               )                           
                          );
                                
                                
                                
                                
                                
                                
            *(im1+count)   = NivGris0+NivGris1+NivGris2+NivGris3;
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
