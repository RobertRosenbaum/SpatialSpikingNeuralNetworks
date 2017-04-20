

/*
 
 Binary search for rheobase of an EIF.
 Searches within the interval I0min and I0max.
 Looks for spike within first Tmax msec.
 Max number of iterations given by niter, doesn't need to be very large.

rheobase=CalcRheoBaseEIF(C,gl,Vl,DeltaT,VT,Vth,v0,I0min,I0max,Tmax,dt,niter);


 */


#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


int i,j,k,m,n,niter,NTmax,flag;
double I0min,I0max,v,Tmax,dt,v0,C,Vleak,DeltaT,VT,gl,Vth,I0,*rheobase;


C=mxGetScalar(prhs[0]);
gl=mxGetScalar(prhs[1]);
Vleak=mxGetScalar(prhs[2]);
DeltaT=mxGetScalar(prhs[3]);
VT=mxGetScalar(prhs[4]);
Vth=mxGetScalar(prhs[5]);

v0=mxGetScalar(prhs[6]);
I0min=mxGetScalar(prhs[7]);
I0max=mxGetScalar(prhs[8]);
Tmax = mxGetScalar(prhs[9]);
dt = mxGetScalar(prhs[10]);
niter=(int)mxGetScalar(prhs[11]);


/******
 * Finished importing variables.
 *******/

NTmax=(int)round(Tmax/dt);


/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
rheobase=mxGetPr(plhs[0]);





v=v0;
    
I0=(I0min+I0max)/2;

for(j=0;j<niter;j++){ 
    
    
    v=v0;
    flag=0;

    for(i=0;i<NTmax && flag==0;i++){
        v+=(I0-gl*(v-Vleak)+gl*DeltaT*exp((v-VT)/DeltaT))*dt/C;
        if(v>=Vth){
            flag=1;
            break;
        }
    }
    if(flag==0){
        I0min=I0;
        I0=(I0+I0max)/2;        
    }
    else{
        I0max=I0;
        I0=(I0+I0min)/2;
    }
    
}


rheobase[0]=I0;


}




