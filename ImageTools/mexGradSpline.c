#include "mex.h"
#include "mat.h"
#include "math.h"

void mexGradSpline(double *Upp,double *Vpp,double *im1x,double *im1y, mwSize m, mwSize n,double *im0, mwSize mc, mwSize nc)
{
    int ii,jj,IX1,IY1,mcnc,ind,count=0;
    double Xi,Yi,Xx,Yy,XxYy,NivGrisx,NivGrisy;
        double a,b,c,d,e,f,g,h,i,j,K,L,M,N,O,p;
double CR[16];
    mcnc=mc*nc;
    for (jj=0;jj<n;jj++){
        for (ii=0;ii<m;ii++){
/*             Xi=*(Upp+count)+ii+1;
             Yi=*(Vpp+count)+jj+1;*/
            Xi=*(Upp+count);
            Yi=*(Vpp+count);
            IX1 = floor(Xi);
            IY1 = floor(Yi);
            Xx=Xi-IX1;
            Yy=Yi-IY1;
            XxYy = Xx * Yy;
            ind=IX1+(IY1-1)*mc-1;
             if ((IX1>1)&&(IY1>1)&&(IX1<mc-3)&&(IY1<nc-3)){
           a = *(im0+ind-1-mc);/*IPixDef(IX1 - 1, IY1 - 1)*/
            b = *(im0+ind-1);/*IPixDef(IX1 - 1, IY1)*/
            c = *(im0+ind-1+mc);/*IPixDef(IX1 - 1, IY1 + 1)*/
            d = *(im0+ind-1+mc+mc);/*IPixDef(IX1 - 1, IY1 + 2)*/
            e = *(im0+ind-mc);/*IPixDef(IX1, IY1 - 1)*/
            f = *(im0+ind);/*IPixDef(IX1, IY1)*/
            g = *(im0+ind+mc);/*IPixDef(IX1, IY1 + 1)*/
            h = *(im0+ind+mc+mc);/*IPixDef(IX1, IY1 + 2)*/
            i = *(im0+ind+1-mc);/*IPixDef(IX1 + 1, IY1 - 1)*/
            j = *(im0+ind+1);/*IPixDef(IX1 + 1, IY1)*/
            K = *(im0+ind+1+mc);/*IPixDef(IX1 + 1, IY1 + 1)*/
            L = *(im0+ind+1+mc+mc);/*IPixDef(IX1 + 1, IY1 + 2)*/
            M = *(im0+ind+2-mc);/*IPixDef(IX1 + 2, IY1 - 1)*/
            N = *(im0+ind+2);/*IPixDef(IX1 + 2, IY1)*/
            O = *(im0+ind+2+mc);/*IPixDef(IX1 + 2, IY1 + 1)*/
            p = *(im0+ind+2+mc+mc);/*IPixDef(IX1 + 2, IY1 + 2)*/

            
            
                       CR[0] = 6.22222222222222E-02 * M + 1.77777777777778E-02 * p -
            0.106666666666667 * O + 2.66666666666667E-02 * N + 0.217777777777778 * a +
            6.22222222222222E-02 * d + 9.33333333333333E-02 * b - 0.373333333333333 * c
            + 0.04 * f - 0.16 * g + 2.66666666666667E-02 * h - 0.16 * j -
            0.106666666666667 * L + 0.64 * K + 9.33333333333333E-02 * e -
            0.373333333333333 * i;
            
            CR[1] = 0.24 * N + 2.66666666666667E-02 * p - 0.16 * O -
            0.373333333333333 * a + 9.33333333333333E-02 * d + 0.84 * b - 0.56 * c +
            0.04 * h + 0.36 * f - 0.24 * g + 0.96 * K - 1.44 * j - 0.16 * L -
            0.106666666666667 * M - 0.16 * e + 0.64 * i;
            
            CR[2] = -0.133333333333333 * N + 0.133333333333333 * O -
            4.44444444444444E-02 * p + 0.155555555555556 * a - 0.466666666666667 * b +
            0.466666666666667 * c - 0.155555555555556 * d - 0.2 * f + 0.2 * g -
            6.66666666666667E-02 * h - 0.8 * K + 0.8 * j + 0.266666666666667 * L +
            4.44444444444444E-02 * M + 6.66666666666667E-02 * e - 0.266666666666667 * i;
            
            CR[3] = -0.133333333333333 * N - 0.2 * f + 0.8 * j -
            0.466666666666667 * b;
            
            CR[4] = -0.16 * L + 0.96 * K + 9.33333333333333E-02 * M +
            2.66666666666667E-02 * p - 0.16 * O + 0.04 * N - 0.373333333333333 * a -
            0.106666666666667 * d - 0.16 * b + 0.64 * c + 0.36 * f - 1.44 * g + 0.24 * h
            - 0.24 * j + 0.84 * e - 0.56 * i;
            
            CR[5] = -0.24 * L - 0.16 * M + 0.36 * N + 0.04 * p - 0.24 * O
            + 0.64 * a - 0.16 * d - 1.44 * b + 0.96 * c + 0.36 * h + 3.24 * f - 2.16 * g
            + 1.44 * K - 2.16 * j - 1.44 * e + 0.96 * i;
            
            CR[6] = 1.2 * j + 0.4 * L - 6.66666666666667E-02 * p +
            6.66666666666667E-02 * M - 0.2 * N + 0.2 * O - 0.266666666666667 * a + 0.8 *
            b - 0.8 * c + 0.266666666666667 * d - 1.8 * f + 1.8 * g - 0.6 * h - 1.2 * K
            + 0.6 * e - 0.4 * i;
            
            CR[7] = -0.2 * N - 1.8 * f + 1.2 * j + 0.8 * b;
            
            CR[8] = -0.155555555555556 * M - 4.44444444444444E-02 * p +
            0.266666666666667 * O - 6.66666666666667E-02 * N + 0.155555555555556 * a +
            4.44444444444444E-02 * d + 6.66666666666667E-02 * b - 0.266666666666667 * c
            - 0.2 * f + 0.8 * g - 0.133333333333333 * h + 0.2 * j + 0.133333333333333 *
            L - 0.8 * K - 0.466666666666667 * e + 0.466666666666667 * i;
            
            CR[9] = 0.266666666666667 * M - 0.6 * N - 6.66666666666667E-02
            * p + 0.4 * O - 0.266666666666667 * a + 6.66666666666667E-02 * d + 0.6 * b -
            0.4 * c - 0.2 * h - 1.8 * f + 1.2 * g - 1.2 * K + 1.8 * j + 0.2 * L + 0.8 *
            e - 0.8 * i;
            
            CR[10] = 0.111111111111111 * p - 0.111111111111111 * M +
            0.333333333333333 * N - 0.333333333333333 * O + 0.111111111111111 * a -
            0.333333333333333 * b + 0.333333333333333 * c - 0.111111111111111 * d + f -
            g + 0.333333333333333 * h + K - j - 0.333333333333333 * L -
            0.333333333333333 * e + 0.333333333333333 * i;
            
            CR[11] = 0.333333333333333 * N - j + f - 0.333333333333333 * b;
            
            CR[12] = -0.466666666666667 * e - 0.2 * f + 0.8 * g -
            0.133333333333333 * h;
            
            CR[13] = 0.8 * e - 0.2 * h - 1.8 * f + 1.2 * g;
            
            CR[14] = -0.333333333333333 * e + f + 0.333333333333333 * h -
            g;
            
            CR[15] =f; 
            

            NivGrisx = Yy * (CR[0] + Yy * (CR[1] + CR[2] * Yy)
            + XxYy * (CR[5] + CR[6] * Yy + CR[9] * Xx + CR[10] * XxYy)
            + Xx * (CR[8] * Xx + CR[4]))
            +XxYy * (
            + Yy * (CR[5] + CR[6] * Yy + CR[9] * Xx + CR[10] * XxYy)
            + XxYy * (CR[9]+ CR[10] * Yy)
            + (CR[8] * Xx + CR[4])
            + Xx * CR[8]
            );
            NivGrisx = NivGrisx + CR[3] + Xx * (2*CR[7] + 3*CR[11] * Xx);
 
            NivGrisy = Xx * (CR[0] + Yy * (CR[1] + CR[2] * Yy)
            + XxYy * (CR[5] + CR[6] * Yy + CR[9] * Xx + CR[10] * XxYy)
            + Xx * (CR[8] * Xx + CR[4]))
            + XxYy * ((CR[1] + CR[2] * Yy)+ Yy * CR[2]
            + Xx * (CR[5] + CR[6] * Yy + CR[9] * Xx + CR[10] * XxYy)
            + XxYy * ( CR[6] +  CR[10] * Xx));
            NivGrisy = NivGrisy + CR[12] + Yy * (2*CR[13] + 3*CR[14] * Yy);
            
            *(im1x+count)   = NivGrisx;
            *(im1y+count)   = NivGrisy;
             }
            else {
             *(im1x+count)   = 0;
            *(im1y+count)   = 0;
               
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

    mexGradSpline(up,vp,im1x,im1y,mrows1,ncols1,im0,mrows0,ncols0);


}
