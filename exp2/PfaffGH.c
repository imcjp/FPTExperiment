#include "mex.h"
#include <math.h>
#define EPS 2.220446049250313e-16
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#define DELTA 1e-5
#define GOLDEN_SECTION_NUMBER 1.618033988749895
#define ROUND(X) (double)((long long)(X+0.5))
double hypergeomBasedPfaff(double x,double y){
    int t=1;
    double a,p,r;
    r=0;
    a=1/(y+1);
    p=1-a;
    while(a>EPS){
        r+=a;
        a*=p*(t/(t+x));
        t++;
    }
    return r;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n,i;
    double *px,*py,*pz,tmp;
    int tp=0;
    if(nrhs < 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin","need two matrice");
    }
    if(nrhs==3){
        tp = * mxGetPr(prhs[2]);
    }
    if((mxGetN(prhs[0]) == mxGetN(prhs[1]))&&(mxGetM(prhs[0]) == mxGetM(prhs[1]))) {
        n=mxGetN(prhs[0])*mxGetM(prhs[0]);
        plhs[0]=mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL);
        px = mxGetPr(prhs[0]);
        py = mxGetPr(prhs[1]);
        pz = mxGetPr(plhs[0]);
        for(i=0;i<n;++i,++px,++py,++pz){
            *pz=hypergeomBasedPfaff(*px,*py);
        }
    }else if((mxGetN(prhs[0]) == 1)||(mxGetM(prhs[0]) == 1)) {
        n=mxGetN(prhs[1])*mxGetM(prhs[1]);
        plhs[0]=mxCreateDoubleMatrix(mxGetM(prhs[1]),mxGetN(prhs[1]),mxREAL);
        px = mxGetPr(prhs[0]);
        tmp=*px;
        py = mxGetPr(prhs[1]);
        pz = mxGetPr(plhs[0]);
        for(i=0;i<n;++i,++py,++pz){
            *pz=hypergeomBasedPfaff(tmp,*py);
        }
    }else if((mxGetN(prhs[1]) ==1)||(mxGetM(prhs[1]) ==1)) {
        n=mxGetN(prhs[0])*mxGetM(prhs[0]);
        plhs[0]=mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL);
        px = mxGetPr(prhs[0]);
        py = mxGetPr(prhs[1]);
        tmp=*py;
        pz = mxGetPr(plhs[0]);
        for(i=0;i<n;++i,++px,++pz){
            *pz=hypergeomBasedPfaff(*px,tmp);
        }
    }else{
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:innerDimensions","Matrix dimensions must agree.");
    }
}