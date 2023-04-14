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

void LSReinit(double *ls,double *lsnew, double dt,double dx, mwSize m, mwSize n ){ 

    int i,j,count;
    double H,Nplus,Nminus;
    double Dxplus,Dyplus,Dzplus;
    double Dxminus,Dyminus,Dzminus;
    double Dxpp,Dxmm,Dxmp,Dxpm;
    double Dypp,Dymm,Dymp,Dypm;
    double Dzpp,Dzmm,Dzmp,Dzpm;
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
            
            if (*(ls+count)<0){
                Dxpp=0.5*Dxplus*(1+SIGN(Dxplus));
                Dxmm=0.5*Dxminus*(1-SIGN(Dxminus));
                
                Dypp=0.5*Dyplus*(1+SIGN(Dyplus));
                Dymm=0.5*Dyminus*(1-SIGN(Dyminus));
                
                Nminus=sqrt(MAX(Dxmm*Dxmm,Dxpp*Dxpp)+MAX(Dymm*Dymm,Dypp*Dypp));
                
                H=-Nminus+1;
                
            }
            else {
                Dxpm=0.5*Dxplus*(1-SIGN(Dxplus));
                Dxmp=0.5*Dxminus*(1+SIGN(Dxminus));
                Dypm=0.5*Dyplus*(1-SIGN(Dyplus));
                Dymp=0.5*Dyminus*(1+SIGN(Dyminus));
                
                
                Nplus=sqrt(MAX(Dxmp*Dxmp,Dxpm*Dxpm)+MAX(Dymp*Dymp,Dypm*Dypm));
                H=Nplus-1;
            }
            
            
            *(lsnew+count)=*(ls+count)-dt*H;
            
            count++;
        }
    }
}
void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *ls,*lsnew,*lstmp;
    mwSize dims[3],ndims;
mxArray *tmp;
double dx,dt,t,tmax;
int i;
    mwSize m,n,p,mn;
    if(nrhs!=3)
        mexErrMsgTxt("Three inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");
    dx=mxGetScalar(prhs[1]);
    dt=0.5*dx;
    dx=1./dx;
    tmax=mxGetScalar(prhs[2]);

    ndims=mxGetNumberOfDimensions(prhs[0]);
    for (i=0;i<ndims;i++){
        dims[i]=MAX(((mxGetDimensions(prhs[0]))[i]),1);
    }
    m=dims[0];
    n=dims[1];
    mn=m*n;
    ls=(double *) mxMalloc(mn*sizeof(double));
    memcpy((void *)(ls),(void *)(mxGetPr(prhs[0])) , mn*sizeof(double));
    
     lsnew=(double *) mxMalloc(mn*sizeof(double));
     
     
     t=0;
     while (t<=tmax){
         LSReinit(ls,lsnew,dt,dx,m,n);
         lstmp=ls;
         ls=lsnew;
         lsnew=lstmp;
      t+=dt;   
     }
 
      plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
   memcpy((void *)(mxGetPr(plhs[0])),(void *)(ls) , mn*sizeof(double));
 
      mxFree(lsnew);
      mxFree(ls);
  
   
     
   
    
}
