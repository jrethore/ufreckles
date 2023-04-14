#include "mex.h"
#include "mat.h"
#include "math.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(SIGN)
#define	SIGN(A)	((A)< 0 ? (-1) : (1))
#endif

void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *ls,*lsnew,*ls0;
    mwSize dims[3],ndims;
    int i,j,k,count;
    double dx,dt,s;
    mwSize m,n,p,mn;
    double H,Nx,Ny,NN;
    double Dx,Dy,Dz;
    if(nrhs!=4)
        mexErrMsgTxt("Four inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");
    dx=mxGetScalar(prhs[2]);
    dx=1./dx;
    dt=mxGetScalar(prhs[3]);
    ls= mxGetPr(prhs[0]);
    ls0= mxGetPr(prhs[1]);
    ndims=mxGetNumberOfDimensions(prhs[0]);
    for (i=0;i<ndims;i++){
        dims[i]=MAX(((mxGetDimensions(prhs[0]))[i]),1);
    }
    m=dims[0];
    n=dims[1];
    mn=m*n;
    plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
    lsnew=mxGetPr(plhs[0]);
    
    count=0;
    
    for (j=0;j<n;j++){

        for (i=0;i<m;i++){
            
            if (i==0){
                Dx=dx * (*(ls+count+1)-*(ls+count));
                Nx=dx*(*(ls0+count+1)-*(ls0+count));
            }
            else if (i==m-1){
                Dx=dx * (*(ls+count)-*(ls+count-1));
                 Nx=dx*(*(ls0+count)-*(ls0+count-1));
           }
            else {
                Dx=0.5*dx * (*(ls+count+1)-*(ls+count-1));
                 Nx=0.5*dx*(*(ls0+count+1)-*(ls0+count-1));
           }
            
            if (j==0){
                Dy=dx * (*(ls+count+m)-*(ls+count));
                 Ny=dx*(*(ls0+count+m)-*(ls0+count));
           }
            else if (j==n-1){
                Dy=dx * (*(ls+count)-*(ls+count-m));
                Ny=dx*(*(ls0+count)-*(ls0+count-m));
           }
            else {
                Dy=0.5*dx * (*(ls+count+m)-*(ls+count-m));
                 Ny=0.5*dx*(*(ls0+count+m)-*(ls0+count-m));
           }
            
            
            
            
            
            NN=1./sqrt(Nx*Nx+Ny*Ny);
            Nx*=NN;
            Ny*=NN;
            s=SIGN(*(ls0+count));
            H=s*(Nx*Dx+Ny*Dy);
            
            
            *(lsnew+count)=*(ls+count)-dt*H;
            
            count++;
        }
    }

    
}
