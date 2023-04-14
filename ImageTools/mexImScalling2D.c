#include "mex.h"
#include "mat.h"
#include "math.h"

void mexScal(double *moy, double *etype, double *im0, mwSize mc, mwSize nc)
{
    int ii,jj,count=0;
    double NivGris;
    double mm,ee;
    mm=0;
    ee=0;
    
    
        for (jj=0;jj<nc;jj++){
            for (ii=0;ii<mc;ii++){
                NivGris=*(im0+count);
                mm+=NivGris;
                count++;
            }
            
        }
   mm*=1./(count);
 
    count=0;
        for (jj=0;jj<nc;jj++){
            for (ii=0;ii<mc;ii++){
                NivGris=*(im0+count)-mm;
                ee+=NivGris*NivGris;
                count++;
            }
            
        }
    ee*=1./(count-1);
    ee=sqrt(ee);
 
    *moy=mm;
    *etype=ee;
}






void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{

    double *im0;
    double *moy,*etype;
    mwSize mrows0,ncols0;

   if(nrhs<1)
        mexErrMsgTxt("Only one input required.");

    if(nlhs!=2)
        mexErrMsgTxt("Two outputs required.");
im0= mxGetData(prhs[0]);
    mrows0 = mxGetM(prhs[0]);
    ncols0 = ((mxGetDimensions(prhs[0]))[1]);
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
moy= mxGetPr(plhs[0]);
etype= mxGetPr(plhs[1]);

    mexScal(moy,etype,im0,mrows0,ncols0);


}
