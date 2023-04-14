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
    double *ls,*lsnew,*v0;
    mwSize dims[3],ndims;
    int i,j,k,count;
    double dx,dt,s;
    mwSize m,n,p,mn;
    double H,Nplus,Nminus;
    double Dxplus,Dyplus,Dzplus;
    double Dxminus,Dyminus,Dzminus;
    double Dxpp,Dxmm,Dxmp,Dxpm;
    double Dypp,Dymm,Dymp,Dypm;
    double Dzpp,Dzmm,Dzmp,Dzpm;
    if(nrhs!=4)
        mexErrMsgTxt("Four inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");
    dx=mxGetScalar(prhs[2]);
    dx=1./dx;
    dt=mxGetScalar(prhs[3]);
    ls= mxGetPr(prhs[0]);
    v0= mxGetPr(prhs[1]);
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
                Dxplus=dx * (*(ls+count+1)-*(ls+count));
                Dxminus=dx * (*(ls+count+1)-*(ls+count));
            }
            else if (i==m-1){
                Dxplus=dx * (*(ls+count)-*(ls+count-1));
                Dxminus=dx * (*(ls+count)-*(ls+count-1));
            }
            else {
                Dxplus=dx * (*(ls+count+1)-*(ls+count));
                Dxminus=dx * (*(ls+count)-*(ls+count-1));
            }
            
            if (j==0){
                Dyplus=dx * (*(ls+count+m)-*(ls+count));
                Dyminus=dx * (*(ls+count+m)-*(ls+count));
            }
            else if (j==n-1){
                Dyplus=dx * (*(ls+count)-*(ls+count-m));
                Dyminus=dx * (*(ls+count)-*(ls+count-m));
            }
            else {
                Dyplus=dx * (*(ls+count+m)-*(ls+count));
                Dyminus=dx * (*(ls+count)-*(ls+count-m));
            }
            
            if (*(v0+count)<0){
                Dxpp=0.5*Dxplus*(1+SIGN(Dxplus));
                Dxmm=0.5*Dxminus*(1-SIGN(Dxminus));
                
                Dypp=0.5*Dyplus*(1+SIGN(Dyplus));
                Dymm=0.5*Dyminus*(1-SIGN(Dyminus));
                
                Nminus=sqrt(MAX(Dxmm*Dxmm,Dxpp*Dxpp)+MAX(Dymm*Dymm,Dypp*Dypp));
                
                H=*(v0+count)*Nminus;
                
            }
            else {
                Dxpm=0.5*Dxplus*(1-SIGN(Dxplus));
                Dxmp=0.5*Dxminus*(1+SIGN(Dxminus));
                Dypm=0.5*Dyplus*(1-SIGN(Dyplus));
                Dymp=0.5*Dyminus*(1+SIGN(Dyminus));
                
                
                Nplus=sqrt(MAX(Dxmp*Dxmp,Dxpm*Dxpm)+MAX(Dymp*Dymp,Dypm*Dypm));
                H=*(v0+count)*Nplus;
            }
 
                        
            
            *(lsnew+count)=*(ls+count)-dt*H;
            
            count++;
        }
    }

    
}
