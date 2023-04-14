#include "mex.h"
#include "mat.h"
#include "math.h"
#include <string.h> /* For memcpy() */
#include <stdio.h>
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */


#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(SIGN)
#define	SIGN(A)	((A)< 0 ? (-1) : (1))
#endif

void LSReinit(double *ls,double *lsnew, double dt,double dx, mwSize m, mwSize n , mwSize o ){
    
    int i,j,k,mn,count;
    double H,Nplus,Nminus;
    double Dxplus,Dyplus,Dzplus;
    double Dxminus,Dyminus,Dzminus;
    double Dxpp,Dxmm,Dxmp,Dxpm;
    double Dypp,Dymm,Dymp,Dypm;
    double Dzpp,Dzmm,Dzmp,Dzpm;
    count=0;
    mn=m*n;
    for (k=0;k<o;k++){
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
                if (k==0){
                    Dzplus=dx * (*(ls+count+mn)-*(ls+count));
                    Dzminus=dx * (*(ls+count+mn)-*(ls+count));
                }
                else if (k==o-1){
                    Dzplus=dx * (*(ls+count)-*(ls+count-mn));
                    Dzminus=dx * (*(ls+count)-*(ls+count-mn));
                }
                else {
                    Dzplus=dx * (*(ls+count+mn)-*(ls+count));
                    Dzminus=dx * (*(ls+count)-*(ls+count-mn));
                }
                
                if (*(ls+count)<0){
                    Dxpp=0.5*Dxplus*(1+SIGN(Dxplus));
                    Dxmm=0.5*Dxminus*(1-SIGN(Dxminus));
                    
                    Dypp=0.5*Dyplus*(1+SIGN(Dyplus));
                    Dymm=0.5*Dyminus*(1-SIGN(Dyminus));
                    
                    Dzpp=0.5*Dzplus*(1+SIGN(Dzplus));
                    Dzmm=0.5*Dzminus*(1-SIGN(Dzminus));
                    
                    Nminus=sqrt(MAX(Dxmm*Dxmm,Dxpp*Dxpp)+MAX(Dymm*Dymm,Dypp*Dypp)+MAX(Dzmm*Dzmm,Dzpp*Dzpp));
                    
                    H=-Nminus+1;
                    
                }
                else {
                    Dxpm=0.5*Dxplus*(1-SIGN(Dxplus));
                    Dxmp=0.5*Dxminus*(1+SIGN(Dxminus));
                    Dypm=0.5*Dyplus*(1-SIGN(Dyplus));
                    Dymp=0.5*Dyminus*(1+SIGN(Dyminus));
                    Dzpm=0.5*Dzplus*(1-SIGN(Dzplus));
                    Dzmp=0.5*Dzminus*(1+SIGN(Dzminus));
                    
                    
                    Nplus=sqrt(MAX(Dxmp*Dxmp,Dxpm*Dxpm)+MAX(Dymp*Dymp,Dypm*Dypm)+MAX(Dzmp*Dzmp,Dzpm*Dzpm));
                    H=Nplus-1;
                }
                
                
                *(lsnew+count)=*(ls+count)-dt*H;
                
                count++;
            }
        }
    }
}
void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *ls,*lsnew,*lstmp;
    mwSize dims[3],ndims;
mxArray *tmp;
double dx,dt;
int i;
    mwSize m,n,p,mn;
    if(nrhs!=3)
        mexErrMsgTxt("Three inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");
    dx=mxGetScalar(prhs[1]);
    dx=1./dx;
    dt=mxGetScalar(prhs[2]);

    ndims=mxGetNumberOfDimensions(prhs[0]);
    for (i=0;i<ndims;i++){
        dims[i]=MAX(((mxGetDimensions(prhs[0]))[i]),1);
    }
    m=dims[0];
    n=dims[1];
    p=dims[2];
    mn=m*n*p;
        ls=mxGetPr(prhs[0]);
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
     lsnew=mxGetPr(plhs[0]);
        LSReinit(ls,lsnew,dt,dx,m,n,p);

}
