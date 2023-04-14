#include "mex.h"
#include "mat.h"
#include "math.h"

void mexPhidf3D(double *phidf,double *N,double *im,int *ipix,int *jpix,int *kpix, mwSize mp, mwSize np, mwSize op, mwSize m, mwSize n, mwSize o,mwSize m1,mwSize n1)
{
    
    int i,j,k,l,mn,pix,count=0;
    double gx,gy,gz;
    mn=m*n;
/*    printf("%d %d %d\n", mp,np,op);
    printf("%d %d %d %d\n", m,n,o,mn);
 **/
            for (k=0;k<op;k++){
            for (j=0;j<np;j++){
            for (i=0;i<mp;i++){
                pix=*(ipix+i)+*(jpix+j)*m+*(kpix+k)*mn; 
 /*               printf("%d %d %d \n", i,j,k);
                printf("%d %d %d %d\n", *(ipix+i),*(jpix+j),*(kpix+k),pix);
  */
                if (i==0) gx=*(im+pix+1)-*(im+pix);
                else if (i==m-1) gx=*(im+pix)-*(im+pix-1);
                else gx=0.5*(*(im+pix+1)-*(im+pix-1));
                if (j==0) gy=*(im+pix+m)-*(im+pix);
                else if (j==n-1) gy=*(im+pix)-*(im+pix-m);
                else gy=0.5*(*(im+pix+m)-*(im+pix-m));
                if (k==0) gz=*(im+pix+mn)-*(im+pix);
                else if (k==o-1) gz=*(im+pix)-*(im+pix-mn);
                else gz=0.5*(*(im+pix+mn)-*(im+pix-mn));

                for (l=0;l<n1;l++){
                   *(phidf+count+l*m1)=*(N+count+l*m1)*gx;
                    *(phidf+count+l*m1+n1*m1)=*(N+count+l*m1)*gy;
                    *(phidf+count+l*m1+2*n1*m1)=*(N+count+l*m1)*gz;
                }
  
                count++;
            }
            }
            }

}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{

    int *ipix,*jpix,*kpix;
    double *im0;
    double *N,*phidf;
    mwSize mrows0,ncols0,pstacks0;
    mwSize mpix,npix,ppix;
    mwSize mrows1,ncols1;
    mwSize dims[2];

   if(nrhs<2)
        mexErrMsgTxt("At least two inputs required.");

    if(nlhs!=1)
        mexErrMsgTxt("Only one output required.");

        im0= mxGetData(prhs[0]);
    N= mxGetData(prhs[1]);
    ipix= mxGetData(prhs[2]);
    jpix= mxGetData(prhs[3]);
    kpix= mxGetData(prhs[4]);
    mpix = mxGetM(prhs[2]);
    npix =mxGetM(prhs[3]);
    ppix =mxGetM(prhs[4]);
    mrows1 = mxGetM(prhs[1]);
    ncols1 = mxGetN(prhs[1]);
     mrows0 = mxGetM(prhs[0]);
    ncols0 = ((mxGetDimensions(prhs[0]))[1]);
    pstacks0 = ((mxGetDimensions(prhs[0]))[2]);
    dims[0]=mrows1;
     dims[1]=3*ncols1;

    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);    
phidf= mxGetData(plhs[0]);

    mexPhidf3D(phidf,N,im0,ipix,jpix,kpix,mpix,npix,ppix,mrows0,ncols0,pstacks0,mrows1,ncols1);


}
