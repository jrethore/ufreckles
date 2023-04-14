#include "mex.h"
#include "mat.h"
#include "math.h"

void mexGradSpline(double *Upp,double *Vpp,double *Wpp,double *gx,double *gy,double *gz, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc, mwSize oc)
{
    int ii,jj,IX1,IY1,IZ1,mcnc,ind,count=0;
    double Xi,Yi,Zi,Xx,Yy,Zz;
    double P0X,P1X,P2X,P3X;
    double P0Y,P1Y,P2Y,P3Y;
    double P0Z,P1Z,P2Z,P3Z;
    double dP0X,dP1X,dP2X,dP3X;
    double dP0Y,dP1Y,dP2Y,dP3Y;
    double dP0Z,dP1Z,dP2Z,dP3Z;
    double dx0,dx1,dx2,dx3;
    double dy0,dy1,dy2,dy3;
    double dz0,dz1,dz2,dz3;
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

            dP0X=0.2*(-7/3+Xx*(4-5*Xx/3))+0.2*Xx*((4-5*Xx/3)+Xx*(-5/3));
            dP1X=0.2*(-1+Xx*(-9+5*Xx))+0.2*Xx*((-9+5*Xx)+Xx*(5));
            dP2X=0.2*(4+Xx*(6-5*Xx))+0.2*Xx*((6-5*Xx)+Xx*(-5));
            dP3X=0.2*(-2/3+Xx*(-1+5/3*Xx))+0.2*Xx*((-1+5/3*Xx)+Xx*(5/3));
            
            P0Y=0.2*Yy*(-7/3+Yy*(4-5*Yy/3));
            P1Y=1+0.2*Yy*(-1+Yy*(-9+5*Yy));
            P2Y=0.2*Yy*(4+Yy*(6-5*Yy));
            P3Y=0.2*Yy*(-2/3+Yy*(-1+5/3*Yy));

            dP0Y=0.2*(-7/3+Yy*(4-5*Yy/3))+0.2*Yy*((4-5*Yy/3)+Yy*(-5/3));
            dP1Y=0.2*(-1+Yy*(-9+5*Yy))+0.2*Yy*((-9+5*Yy)+Yy*(5));
            dP2Y=0.2*(4+Yy*(6-5*Yy))+0.2*Yy*((6-5*Yy)+Yy*(-5));
            dP3Y=0.2*(-2/3+Yy*(-1+5/3*Yy))+0.2*Yy*((-1+5/3*Yy)+Yy*(5/3));
            
            P0Z=0.2*Zz*(-7/3+Zz*(4-5*Zz/3));
            P1Z=1+0.2*Zz*(-1+Zz*(-9+5*Zz));
            P2Z=0.2*Zz*(4+Zz*(6-5*Zz));
            P3Z=0.2*Zz*(-2/3+Zz*(-1+5/3*Zz));
            
            dP0Z=0.2*(-7/3+Zz*(4-5*Zz/3))+0.2*Zz*((4-5*Zz/3)+Zz*(-5/3));
            dP1Z=0.2*(-1+Zz*(-9+5*Zz))+0.2*Zz*((-9+5*Zz)+Zz*(5));
            dP2Z=0.2*(4+Zz*(6-5*Zz))+0.2*Zz*((6-5*Zz)+Zz*(-5));
            dP3Z=0.2*(-2/3+Zz*(-1+5/3*Zz))+0.2*Zz*((-1+5/3*Zz)+Zz*(5/3));
                      
            
            dx0=dP0X*(
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
                      
            dx1=dP1X*(
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
                                
            dx2=dP2X*(
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
                                
            dx3=dP3X*(
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
                                
                                
                                
                                
                                
                                
            *(gx+count)   = dx0+dx1+dx2+dx3;
            
            dy0=P0X*(
            dP0Y*(
            P0Z*(*(im0+ind-1+ijmkm))+P1Z*(*(im0+ind-1+ijmk))+P2Z*(*(im0+ind-1+ijmkp))+P3Z*(*(im0+ind-1+ijmkpp))
            )
            +dP1Y*(
            P0Z*(*(im0+ind-1+ijkm))+P1Z*(*(im0+ind-1+ijk))+P2Z*(*(im0+ind-1+ijkp))+P3Z*(*(im0+ind-1+ijkpp))
            )
            +dP2Y*(
            P0Z*(*(im0+ind-1+ijpkm))+P1Z*(*(im0+ind-1+ijpk))+P2Z*(*(im0+ind-1+ijpkp))+P3Z*(*(im0+ind-1+ijpkpp))
            )
            +dP3Y*(
            P0Z*(*(im0+ind-1+ijppkm))+P1Z*(*(im0+ind-1+ijppk))+P2Z*(*(im0+ind-1+ijppkp))+P3Z*(*(im0+ind-1+ijppkpp))
            )
            );
            
            dy1=P1X*(
            dP0Y*(
            P0Z*(*(im0+ind+ijmkm))+P1Z*(*(im0+ind+ijmk))+P2Z*(*(im0+ind+ijmkp))+P3Z*(*(im0+ind+ijmkpp))
            )
            +dP1Y*(
            P0Z*(*(im0+ind+ijkm))+P1Z*(*(im0+ind+ijk))+P2Z*(*(im0+ind+ijkp))+P3Z*(*(im0+ind+ijkpp))
            )
            +dP2Y*(
            P0Z*(*(im0+ind+ijpkm))+P1Z*(*(im0+ind+ijpk))+P2Z*(*(im0+ind+ijpkp))+P3Z*(*(im0+ind+ijpkpp))
            )
            +dP3Y*(
            P0Z*(*(im0+ind+ijppkm))+P1Z*(*(im0+ind+ijppk))+P2Z*(*(im0+ind+ijppkp))+P3Z*(*(im0+ind+ijppkpp))
            )
            );
            
            dy2=P2X*(
            dP0Y*(
            P0Z*(*(im0+ind+1+ijmkm))+P1Z*(*(im0+ind+1+ijmk))+P2Z*(*(im0+ind+1+ijmkp))+P3Z*(*(im0+ind+1+ijmkpp))
            )
            +dP1Y*(
            P0Z*(*(im0+ind+1+ijkm))+P1Z*(*(im0+ind+1+ijk))+P2Z*(*(im0+ind+1+ijkp))+P3Z*(*(im0+ind+1+ijkpp))
            )
            +dP2Y*(
            P0Z*(*(im0+ind+1+ijpkm))+P1Z*(*(im0+ind+1+ijpk))+P2Z*(*(im0+ind+1+ijpkp))+P3Z*(*(im0+ind+1+ijpkpp))
            )
            +dP3Y*(
            P0Z*(*(im0+ind+1+ijppkm))+P1Z*(*(im0+ind+1+ijppk))+P2Z*(*(im0+ind+1+ijppkp))+P3Z*(*(im0+ind+1+ijppkpp))
            )
            );
            
            dy3=P3X*(
            dP0Y*(
            P0Z*(*(im0+ind+2+ijmkm))+P1Z*(*(im0+ind+2+ijmk))+P2Z*(*(im0+ind+2+ijmkp))+P3Z*(*(im0+ind+2+ijmkpp))
            )
            +dP1Y*(
            P0Z*(*(im0+ind+2+ijkm))+P1Z*(*(im0+ind+2+ijk))+P2Z*(*(im0+ind+2+ijkp))+P3Z*(*(im0+ind+2+ijkpp))
            )
            +dP2Y*(
            P0Z*(*(im0+ind+2+ijpkm))+P1Z*(*(im0+ind+2+ijpk))+P2Z*(*(im0+ind+2+ijpkp))+P3Z*(*(im0+ind+2+ijpkpp))
            )
            +dP3Y*(
            P0Z*(*(im0+ind+2+ijppkm))+P1Z*(*(im0+ind+2+ijppk))+P2Z*(*(im0+ind+2+ijppkp))+P3Z*(*(im0+ind+2+ijppkpp))
            )
            );
            
            *(gy+count)   = dy0+dy1+dy2+dy3;
            
                        dz0=P0X*(
                           P0Y*(
                                dP0Z*(*(im0+ind-1+ijmkm))+dP1Z*(*(im0+ind-1+ijmk))+dP2Z*(*(im0+ind-1+ijmkp))+dP3Z*(*(im0+ind-1+ijmkpp))
                               )
                           +P1Y*(
                                dP0Z*(*(im0+ind-1+ijkm))+dP1Z*(*(im0+ind-1+ijk))+dP2Z*(*(im0+ind-1+ijkp))+dP3Z*(*(im0+ind-1+ijkpp))
                               )
                           +P2Y*(
                                dP0Z*(*(im0+ind-1+ijpkm))+dP1Z*(*(im0+ind-1+ijpk))+dP2Z*(*(im0+ind-1+ijpkp))+dP3Z*(*(im0+ind-1+ijpkpp))
                               )
                           +P3Y*(
                                dP0Z*(*(im0+ind-1+ijppkm))+dP1Z*(*(im0+ind-1+ijppk))+dP2Z*(*(im0+ind-1+ijppkp))+dP3Z*(*(im0+ind-1+ijppkpp))
                               )                           
                          );
                      
            dz1=P1X*(
                           P0Y*(
                                dP0Z*(*(im0+ind+ijmkm))+dP1Z*(*(im0+ind+ijmk))+dP2Z*(*(im0+ind+ijmkp))+dP3Z*(*(im0+ind+ijmkpp))
                               )
                           +P1Y*(
                                dP0Z*(*(im0+ind+ijkm))+dP1Z*(*(im0+ind+ijk))+dP2Z*(*(im0+ind+ijkp))+dP3Z*(*(im0+ind+ijkpp))
                               )
                           +P2Y*(
                                dP0Z*(*(im0+ind+ijpkm))+dP1Z*(*(im0+ind+ijpk))+dP2Z*(*(im0+ind+ijpkp))+dP3Z*(*(im0+ind+ijpkpp))
                               )
                           +P3Y*(
                                dP0Z*(*(im0+ind+ijppkm))+dP1Z*(*(im0+ind+ijppk))+dP2Z*(*(im0+ind+ijppkp))+dP3Z*(*(im0+ind+ijppkpp))
                               )                           
                          );
                                
            dz2=P2X*(
                           P0Y*(
                                dP0Z*(*(im0+ind+1+ijmkm))+dP1Z*(*(im0+ind+1+ijmk))+dP2Z*(*(im0+ind+1+ijmkp))+dP3Z*(*(im0+ind+1+ijmkpp))
                               )
                           +P1Y*(
                                dP0Z*(*(im0+ind+1+ijkm))+dP1Z*(*(im0+ind+1+ijk))+dP2Z*(*(im0+ind+1+ijkp))+dP3Z*(*(im0+ind+1+ijkpp))
                               )
                           +P2Y*(
                                dP0Z*(*(im0+ind+1+ijpkm))+dP1Z*(*(im0+ind+1+ijpk))+dP2Z*(*(im0+ind+1+ijpkp))+dP3Z*(*(im0+ind+1+ijpkpp))
                               )
                           +P3Y*(
                                dP0Z*(*(im0+ind+1+ijppkm))+dP1Z*(*(im0+ind+1+ijppk))+dP2Z*(*(im0+ind+1+ijppkp))+dP3Z*(*(im0+ind+1+ijppkpp))
                               )                           
                          );
                                
            dz3=P3X*(
                           P0Y*(
                                P0Z*(*(im0+ind+2+ijmkm))+P1Z*(*(im0+ind+2+ijmk))+dP2Z*(*(im0+ind+2+ijmkp))+dP3Z*(*(im0+ind+2+ijmkpp))
                               )
                           +P1Y*(
                                P0Z*(*(im0+ind+2+ijkm))+P1Z*(*(im0+ind+2+ijk))+dP2Z*(*(im0+ind+2+ijkp))+dP3Z*(*(im0+ind+2+ijkpp))
                               )
                           +P2Y*(
                                P0Z*(*(im0+ind+2+ijpkm))+P1Z*(*(im0+ind+2+ijpk))+dP2Z*(*(im0+ind+2+ijpkp))+dP3Z*(*(im0+ind+2+ijpkpp))
                               )
                           +P3Y*(
                                P0Z*(*(im0+ind+2+ijppkm))+P1Z*(*(im0+ind+2+ijppk))+dP2Z*(*(im0+ind+2+ijppkp))+dP3Z*(*(im0+ind+2+ijppkpp))
                               )                           
                          );

            
            *(gz+count)   = dz0+dz1+dz2+dz3;
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
