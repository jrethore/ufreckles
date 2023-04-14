#include "mex.h"
#include "mat.h"
#include "math.h"

void mexInv66(double *q,double *iq, mwSize mc, mwSize nc)
{
    int m,n1,ll,ii,jj,ind,count=0;
    int k,is,it,ic;
    int n,nn,n2;
    double W,B,C,AINV[6][6],A[6][12],ID[6];
    n=6;
    nn=n;
    n2=2*n;
    for (ll=0;ll<nc;ll++){
        ind=0;

        for (jj=0;jj<6;jj++){
            for (ii=0;ii<6;ii++){
                A[ii][jj]=*(q+count+ind);
                ind++;
            }
        }
        
        for (ii=0;ii<6;ii++){
            for (jj=n;jj<n2;jj++){
                A[ii][jj]=0.;
            }
        }
        for (ii=0;ii<6;ii++){
            A[ii][n+ii]=1.;
            ID[ii]=ii;
        }
        
        for (k=0;k<n-1;k++){
            
            is=k;
            it=k;
            B=fabs(A[k][k]);
            for (ii=k;ii<n;ii++){
                for (jj=k;jj<n;jj++){
                    if (fabs(A[ii][jj])>B){
                        is=ii;
                        it=jj;
                        B=fabs(A[ii][jj]);
                    }
                }
            }
            if (is>k){
                for (jj=k;jj<n2;jj++){
                    C=A[is][jj];
                    A[is][jj]=A[k][jj];
                    A[k][jj]=C;
                }
            }
            if (it>k){
                ic=ID[k];
                ID[k]=ID[it];
                ID[it]=ic;
                for (ii=0;ii<n;ii++){
                    C=A[ii][it];
                    A[ii][it]=A[ii][k];
                    A[ii][k]=C;
                }
            }
            if (A[k][k]==0) {
                mexPrintf("Pivot %d %d %g\n",k+1,k+1,A[k][k]);
                ind=0;
                for (ii=0;ii<6;ii++){
                    for (jj=0;jj<6;jj++){
                        mexPrintf("     %7.4f",*(q+count+ind));
                        ind++;
                    }
                    mexPrintf("\n");
                }
                
                mexErrMsgTxt("Singular Matrix.");
            }
            for (jj=k+1;jj<n2;jj++){
                A[k][jj]*=(1./A[k][k]);
                for (ii=k+1;ii<n;ii++){
                    W=A[ii][k]*A[k][jj];
                    A[ii][jj]+=-W;
                    if (fabs(A[ii][jj])<0.000000001*fabs(W)){
                        A[ii][jj]=0.;
                    }
                }
            }

        }
        
        if (A[n-1][n-1]==0) {
            mexPrintf("Pivot %d %d %g\n",n,n,A[n-1][n-1]);
            ind=0;
            for (ii=0;ii<6;ii++){
                for (jj=0;jj<6;jj++){
                    mexPrintf("     %7.4f",*(q+count+ind));
                    ind++;
                }
                mexPrintf("\n");
            }
            mexErrMsgTxt("Singular Matrix.");
        }
        
        for (jj=n;jj<n2;jj++){
            A[n-1][jj]*=1./A[n-1][n-1];
        }
        for (m=1;m<n;m++){
            ii=n-1-m;
            for (k=n-m;k<n;k++){
                for (jj=n;jj<n2;jj++){
                    A[ii][jj]+=-A[ii][k]*A[k][jj];
                }
            }
        }
        for (ii=0;ii<6;ii++){
            for (jj=0;jj<6;jj++){
                if (ID[jj]==ii){
                    for (k=n;k<n2;k++){
                        AINV[ii][k-n]=A[jj][k];
                    }
                }
            }
        }
        
        
        ind=0;
        for (jj=0;jj<6;jj++){
            for (ii=0;ii<6;ii++){
                *(iq+count+ind)=AINV[ii][jj];
                ind++;
                
            }
        }
        count+=36;
    }
    
}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    
    double *q,*iq;
    mwSize mrows0,ncols0;
    mwSize dims[2];
    
    
    
    q= mxGetData(prhs[0]);
    mrows0 = mxGetM(prhs[0]);
    ncols0 = mxGetN(prhs[0]);
    dims[0]=mrows0;dims[1]=ncols0;;
    if(mrows0!=36)
        mexErrMsgTxt("Only 6x6 matrix !");
    if(nrhs>1)
        mexErrMsgTxt("Only one input required.");
    
    if(nlhs!=1)
        mexErrMsgTxt("Only one output required.");
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    iq= mxGetData(plhs[0]);
    
    mexInv66(q,iq,mrows0,ncols0);
    
    
}
