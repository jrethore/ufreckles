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
    double H,Nx,Ny,Nz,NN;
    double Dxplus,Dyplus,Dzplus;
    double Dxminus,Dyminus,Dzminus;
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
    p=dims[2];
    mn=m*n;
     plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
    lsnew=mxGetPr(plhs[0]);
    
    count=0;
        for (k=0;k<p;k++){

    for (j=0;j<n;j++){

        for (i=0;i<m;i++){
            
            if (i==0){
                Dxplus=dx * (*(ls+count+1)-*(ls+count));
                Dxminus=dx * (*(ls+count+1)-*(ls+count));
                Nx=dx*(*(ls0+count+1)-*(ls0+count));
            }
            else if (i==m-1){
                Dxplus=dx * (*(ls+count)-*(ls+count-1));
                Dxminus=dx * (*(ls+count)-*(ls+count-1));
                 Nx=dx*(*(ls0+count)-*(ls0+count-1));
           }
            else {
                Dxplus=dx * (*(ls+count+1)-*(ls+count));
                Dxminus=dx * (*(ls+count)-*(ls+count-1));
                 Nx=0.5*dx*(*(ls0+count+1)-*(ls0+count-1));
           }
            
            if (j==0){
                Dyplus=dx * (*(ls+count+m)-*(ls+count));
                Dyminus=dx * (*(ls+count+m)-*(ls+count));
                 Ny=dx*(*(ls0+count+m)-*(ls0+count));
           }
            else if (j==n-1){
                Dyplus=dx * (*(ls+count)-*(ls+count-m));
                Dyminus=dx * (*(ls+count)-*(ls+count-m));
                Ny=dx*(*(ls0+count)-*(ls0+count-m));
           }
            else {
                Dyplus=dx * (*(ls+count+m)-*(ls+count));
                Dyminus=dx * (*(ls+count)-*(ls+count-m));
                 Ny=0.5*dx*(*(ls0+count+m)-*(ls0+count-m));
           }
            
            if (k==0){
                Dzplus=dx * (*(ls+count+mn)-*(ls+count));
                Dzminus=dx * (*(ls+count+mn)-*(ls+count));
                 Nz=dx*(*(ls0+count+mn)-*(ls0+count));
           }
            else if (k==p-1){
                Dzplus=dx * (*(ls+count)-*(ls+count-mn));
                Dzminus=dx * (*(ls+count)-*(ls+count-mn));
                Nz=dx*(*(ls0+count)-*(ls0+count-mn));
           }
            else {
                Dzplus=dx * (*(ls+count+mn)-*(ls+count));
                Dzminus=dx * (*(ls+count)-*(ls+count-mn));
                 Nz=0.5*dx*(*(ls0+count+mn)-*(ls0+count-mn));
           }
            
            
            
            
            
            NN=1./sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
            Nx*=NN;
            Ny*=NN;
            Nz*=NN;
            s=SIGN(*(ls0+count));
            H=MAX(s*Nx,0)*Dxminus+MIN(s*Nx,0)*Dxplus
            +MAX(s*Ny,0)*Dyminus+MIN(s*Ny,0)*Dyplus
            +MAX(s*Nz,0)*Dzminus+MIN(s*Nz,0)*Dzplus;
            
            
            *(lsnew+count)=*(ls+count)-dt*H;
            
            count++;
        }
    }
        }

    
}
