

/* To use function from matlab, first compile by entering this into the Matlab command window:
   mex LIFNetworkTimeDep0.c
   Then call the function like this:
 [s,Ix,Ie,Ii,v0]=EIFTwoPopNetwork(Ix1e,Ix2e,Ix1i,Ix2i,Ne,Ni,Ne1,Ni1,Jex,Jix,Jee,Jei,Jie,Jii,rxe,rxi,Kee,Kei,Kie,Kii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynx,tausyne,tausyni,V0,T,dt,maxns,Irecord);

 Ix1e, Ix2e, Ix1i, Ix2i are inputs to excitatory and inhibitory (e/i) neurons
 in populations 1 and 2.
 
 See EIF2DSpatialNetwork.c for description of other variables
 
 */


#include "mex.h"
#include "math.h"
#include "matrix.h"


#include "mex.h"
#include "math.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifndef srand48
#   define srand48(s) srand(s)
#endif

#ifndef drand48
#   define drand48() (((double)rand())/((double)RAND_MAX))
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


int N,Ne,Ni,Ne1,Ni1,i,j,k,jj,m1,m2,Nt,Kee,Kei,Kie,Kii,maxns,Nrecord,Ntref[2],ns,flag;
double dt,*s,*v,*v0,*JnextE,*JnextI,*alphax,*alphae,*alphai,tausyne,tausyni,*alphaxr,*alphaer,*alphair,*vr,*Ix1e,*Ix2e,*Ix1i,*Ix2i,Jex,Jix,tausynx,rxe,rxi;
double Jee,Jei,Jie,Jii,T,*Irecord,*C,*Vleak,*DeltaT,*VT,*tref,*gl,*Vth,*Vre,*Vlb;
int *Wee,*Wei,*Wie,*Wii,*refstate;

/*t
    uint64_t        start,start1,start2,start3,start4,start5;
    uint64_t        end,end1,end2,end3,end4,end5;
    uint64_t        elapsed,elapsed1,elapsed2,elapsed3,elapsed4,elapsed5;

    elapsed=0; elapsed1=0; elapsed2=0; elapsed3=0; elapsed4=0; elapsed5=0;
t*/   


/******
 * Import variables from matlab
 * This is messy looking and is specific to mex.
 * Ignore if you're implementing this outside of mex.
 *******/
Ix1e = mxGetPr(prhs[0]);
Nt = mxGetM(prhs[0]);
m1 = mxGetN(prhs[0]);
if(Nt==1 && m1!=1){
    Nt=m1;
    m1=1;
}
if(Nt==1 || m1!=1)
    mexErrMsgTxt("Ix1e should be Ntx1 or 1xNt");

Ix2e = mxGetPr(prhs[1]);
m1 = mxGetM(prhs[1]);
m2 = mxGetN(prhs[1]);
if(!((m1==Nt && m2==1)||(m1==1 && m2==Nt))){
   mexPrintf("\n%d %d %d\n",m1,m2,Nt);    
   mexErrMsgTxt("Ix2e should be Ntx1 or 1xNt");
}

Ix1i = mxGetPr(prhs[2]);
m1 = mxGetM(prhs[2]);
m2 = mxGetN(prhs[2]);
if(!((m1==Nt && m2==1)||(m1==1 && m2==Nt)))
   mexErrMsgTxt("Ix1i should be Ntx1 or 1xNt");

Ix2i = mxGetPr(prhs[3]);
m1 = mxGetM(prhs[3]);
m2 = mxGetN(prhs[3]);
if(!((m1==Nt && m2==1)||(m1==1 && m2==Nt)))
   mexErrMsgTxt("Ix2i should be Ntx1 or 1xNt");


Ne=(int)mxGetScalar(prhs[4]);
Ni=(int)mxGetScalar(prhs[5]);

Ne1=(int)mxGetScalar(prhs[6]);
Ni1=(int)mxGetScalar(prhs[7]);

Jex=mxGetScalar(prhs[8]);
Jix=mxGetScalar(prhs[9]);

Jee= mxGetScalar(prhs[10]);
Jei= mxGetScalar(prhs[11]);
Jie= mxGetScalar(prhs[12]);
Jii= mxGetScalar(prhs[13]);

rxe=(int)mxGetScalar(prhs[14]);
rxi=(int)mxGetScalar(prhs[15]);

Kee=(int)mxGetScalar(prhs[16]);
Kei=(int)mxGetScalar(prhs[17]);
Kie=(int)mxGetScalar(prhs[18]);
Kii=(int)mxGetScalar(prhs[19]);



C=mxGetPr(prhs[20]);
m1 = mxGetN(prhs[20]);
m2 = mxGetM(prhs[20]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
gl=mxGetPr(prhs[21]);
m1 = mxGetN(prhs[21]);
m2 = mxGetM(prhs[21]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vleak=mxGetPr(prhs[22]);
m1 = mxGetN(prhs[22]);
m2 = mxGetM(prhs[22]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
DeltaT=mxGetPr(prhs[23]);
m1 = mxGetN(prhs[23]);
m2 = mxGetM(prhs[23]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
VT=mxGetPr(prhs[24]);
m1 = mxGetN(prhs[24]);
m2 = mxGetM(prhs[24]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
tref=mxGetPr(prhs[25]);
m1 = mxGetN(prhs[25]);
m2 = mxGetM(prhs[25]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vth=mxGetPr(prhs[26]);
m1 = mxGetN(prhs[26]);
m2 = mxGetM(prhs[26]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vre=mxGetPr(prhs[27]);
m1 = mxGetN(prhs[27]);
m2 = mxGetM(prhs[27]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vlb=mxGetPr(prhs[28]);
m1 = mxGetN(prhs[28]);
m2 = mxGetM(prhs[28]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");



tausynx=mxGetScalar(prhs[29]);
tausyne=mxGetScalar(prhs[30]);
tausyni=mxGetScalar(prhs[31]);

v0 = mxGetPr(prhs[32]);
N = mxGetM(prhs[32]);
m2 = mxGetN(prhs[32]);
if(N==1 && m2!=1)
    N=m2;

T = mxGetScalar(prhs[33]);
dt = mxGetScalar(prhs[34]);

maxns = ((int)mxGetScalar(prhs[35]));

Irecord=mxGetPr(prhs[36]);
Nrecord = mxGetN(prhs[36]);
m1 = mxGetM(prhs[36]);
if(Nrecord==1 && m1!=1){
    Nrecord=m1;
    m1=1;
}
if(Nrecord==1 || m1!=1)
    mexErrMsgTxt("Irecord should be Nrecordx1 or 1xNrecord");


/******
 * Finished importing variables.
 *******/


N=Ne+Ni;

Nt=(int)round(T/dt);

/******
 * Now allocate new variables.
 * This is also mex specific.  Use malloc in C, etc.
 *****/

/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(2, maxns, mxREAL);
s=mxGetPr(plhs[0]);

/* Allocate output vector */
plhs[1] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
alphaxr=mxGetPr(plhs[1]);

/* Allocate output vector */
plhs[2] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
alphaer=mxGetPr(plhs[2]);

/* Allocate output vector */
plhs[3] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
alphair=mxGetPr(plhs[3]);

/* Allocate output vector */
plhs[4] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
vr=mxGetPr(plhs[4]);


/* Allocate membrane potential */
v=mxMalloc(N*sizeof(double));
JnextE=mxMalloc(N*sizeof(double));
JnextI=mxMalloc(N*sizeof(double));
alphae=mxMalloc(N*sizeof(double));
alphai=mxMalloc(N*sizeof(double));
alphax=mxMalloc(N*sizeof(double));
Wee=mxMalloc(Ne*Kee*sizeof(int));
Wei=mxMalloc(Ni*Kei*sizeof(int));
Wie=mxMalloc(Ne*Kie*sizeof(int));
Wii=mxMalloc(Ni*Kii*sizeof(int));
refstate=mxMalloc(N*sizeof(int));

/*****
 * Finished allocating variables
 ****/

/* Inititalize v */
for(j=0;j<N;j++){
    v[j]=v0[j]; 
    refstate[j]=0;
    JnextE[j]=0;
    JnextI[j]=0;
    alphae[j]=0;
    alphai[j]=0;
    alphax[j]=0;
   }


for(jj=0;jj<Nrecord;jj++){
      /* Find index into local variables */
    
        /* Find index into local variables */
        j=(int)round(Irecord[jj]-1);
                
      if(j>=N || j<0)
         mexErrMsgTxt("Irecord contains out of bounds indices.");
  
      alphaer[jj+Nrecord*0]=alphae[j];
      alphair[jj+Nrecord*0]=alphai[j];
      alphaxr[jj+Nrecord*0]=alphax[j];
      vr[jj+Nrecord*0]=v[j];
}
    

/*t
start = mach_absolute_time();
t*/
/* Initialize connection matrix */
for(j=0;j<Ne;j++){                
    for(k=0;k<Kee;k++)
        Wee[j*Kee+k]=(int)floor(drand48()*((double)(Ne)));
    for(k=0;k<Kie;k++)
        Wie[j*Kie+k]=Ne+(int)floor(drand48()*((double)(Ni)));         
     }

for(j=0;j<Ni;j++){                         
    for(k=0;k<Kei;k++)
        Wei[j*Kei+k]=(int)floor(drand48()*((double)(Ne)));
    for(k=0;k<Kii;k++)
        Wii[j*Kii+k]=Ne+(int)floor(drand48()*((double)(Ni)));         
     }


/*t
end = mach_absolute_time();
elapsed = end - start;
mexPrintf("\nBuilding Ws: %" PRIu64 "\n", elapsed);
t*/



Ntref[0]=(int)round(tref[0]/dt);
Ntref[1]=(int)round(tref[1]/dt);

/* Number of x spikes in whole network at each time bin */
/*
nex=(int)round(((double)Ne)*rxe*dt);
nix=(int)round(((double)Ni)*rxi*dt);
if(nex<100 || nix<100)
    mexWarnMsgTxt("nex or nix is small.  Bad approx to Poisson for input spikes.");
*/

/* Initialize number of spikes */
ns=0;


flag=0;


/*t start = mach_absolute_time(); t*/
/* Time loop */
/* Exit loop and issue a warning if max number of spikes is exceeded */
for(i=1;i<Nt && ns<maxns;i++){
    
    
    /*t start1 = mach_absolute_time(); t*/
    /*t end1 = mach_absolute_time();
     elapsed1 += end1 - start1; t*/   
    /*t start5 = mach_absolute_time();t*/
     for(j=0;j<N;j++){    
    
         alphae[j]-=alphae[j]*(dt/tausyne);
         alphai[j]-=alphai[j]*(dt/tausyni);  
         alphax[j]-=alphax[j]*(dt/tausynx);  
      
         if(j<Ne){
             if(refstate[j]<=0){
                v[j]+=(alphae[j]+alphai[j]+alphax[j]-gl[0]*(v[j]-Vleak[0])+gl[0]*DeltaT[0]*exp((v[j]-VT[0])/DeltaT[0]))*dt/C[0];
                if(j<Ne1)
                    v[j]+=Ix1e[i]*dt/C[0];
                else
                    v[j]+=Ix2e[i]*dt/C[0];                
                v[j]=fmax(v[j],Vlb[0]);
             }
             else{                 
                if(refstate[j]>1)
                   v[j]=Vth[0];/*-=(Vth[0]-Vre[0])/((double)Ntref[0]);*/
                else
                   v[j]=Vre[0];
                refstate[j]--;
             }
             /* If a spike occurs */
             if(v[j]>=Vth[0] && refstate[j]<=0 && ns<maxns){

                 /*t start2 = mach_absolute_time();t*/
                 
                  
                  refstate[j]=Ntref[0];
                  v[j]=Vth[0];       /* reset membrane potential */
                  s[0+2*ns]=i*dt; /* spike time */
                  s[1+2*ns]=j+1;     /* neuron index 1 */
                  ns++;           /* update total number of spikes */


                  /* For each postsynaptic target */
                  for(k=0;k<Kee;k++)
                         JnextE[Wee[j*Kee+k]]+=Jee;
                  for(k=0;k<Kie;k++)
                         JnextE[Wie[j*Kie+k]]+=Jie;

                 /*t end2 = mach_absolute_time();
                  elapsed2 += end2 - start2; t*/
              }                 
         }         
         else{ /* If cell is inhibitory */
            
             if(refstate[j]<=0){
                v[j]+=(alphae[j]+alphai[j]+alphax[j]-gl[1]*(v[j]-Vleak[1])+gl[1]*DeltaT[1]*exp((v[j]-VT[1])/DeltaT[1]))*dt/C[1];
                if(j<Ne+Ni1)
                    v[j]+=Ix1i[i]*dt/C[1];
                else
                    v[j]+=Ix2i[i]*dt/C[1];                
                v[j]=fmax(v[j],Vlb[1]);
             }
             else{                 
                if(refstate[j]>1)
                   v[j]=Vth[1];/*-=(Vth[1]-Vre[1])/((double)Ntref[1]);*/
                else
                   v[j]=Vre[1];
                refstate[j]--;
             }
             
             
              /* If a spike occurs */
              if(v[j]>=Vth[1] && refstate[j]<=0 && ns<maxns){                                                            
                  
                 /*t start2 = mach_absolute_time(); t*/
                  
                  refstate[j]=Ntref[1];
                  v[j]=Vth[1];       /* reset membrane potential */
                  s[0+2*ns]=i*dt; /* spike time */
                  s[1+2*ns]=j+1;     /* neuron index 1 */
                  
                  ns++;           /* update total number of spikes */


                 /* For each postsynaptic target */
                for(k=0;k<Kei;k++)
                        JnextI[Wei[(j-Ne)*Kei+k]]+=Jei;
                 for(k=0;k<Kii;k++)
                        JnextI[Wii[(j-Ne)*Kii+k]]+=Jii;
                
                  /*t end2 = mach_absolute_time();
                  elapsed2 += end2 - start2; t*/
               }
           
             
           }
          

         
        }
     /*t end5 = mach_absolute_time();
     elapsed5 += end5 - start5;       t*/
  

    
    
    
           /* Store recorded variables */
        /*t start3 = mach_absolute_time();   t*/
        for(jj=0;jj<Nrecord;jj++){
            
            /* Find index into local variables */
            j=(int)round(Irecord[jj]-1);

            if(j>=N || j<0)
                mexErrMsgTxt("Bad index in Irecord.");
            
              alphaer[jj+Nrecord*i]=alphae[j];
              alphair[jj+Nrecord*i]=alphai[j];
              alphaxr[jj+Nrecord*i]=alphax[j];
              vr[jj+Nrecord*i]=v[j];
              
              /*t end3 = mach_absolute_time();
              elapsed3 += end3 - start3; t*/               
        }    

        /* Propagate spikes */
       /*t start4 = mach_absolute_time();           t*/
        for(j=0;j<N;j++){                     
          alphae[j]+=JnextE[j]/tausyne;
          alphai[j]+=JnextI[j]/tausyni;
          JnextE[j]=0;
          JnextI[j]=0;            
        }    
     
        for(j=0;j<Ne;j++)
            if(drand48()<rxe*dt)
                alphax[j]+=Jex/tausynx;
     
        for(j=Ne;j<N;j++)
            if(drand48()<rxi*dt)
                alphax[j]+=Jix/tausynx;

     
       /*t end4 = mach_absolute_time();
        elapsed4 += end4 - start4;     t*/
    
}


/*t
end = mach_absolute_time();
elapsed = end - start;
mexPrintf("\nEntire time loop: %" PRIu64 "\n", elapsed);
mexPrintf("\nExternal spikes: %" PRIu64 "\n", elapsed1);
mexPrintf("\nRecurrent spikes: %" PRIu64 "\n", elapsed2);
mexPrintf("\nRecording stuff: %" PRIu64 "\n", elapsed3);
mexPrintf("\nUpdating Jnexts: %" PRIu64 "\n", elapsed4);
mexPrintf("\nj Loops: %" PRIu64 "\n", elapsed5);
mexPrintf("\nj Loops without spike prop: %" PRIu64 "\n", elapsed5-elapsed2);


mexPrintf("\nEntire time loop: %f\n", ((double)(elapsed))/1000000000);
mexPrintf("\nExternal spikes: %f\n", ((double)(elapsed1))/1000000000);
mexPrintf("\nRecurrent spikes: %f\n", ((double)(elapsed2))/1000000000);
mexPrintf("\nRecording stuff: %f\n", ((double)(elapsed3))/1000000000);
mexPrintf("\nUpdating Jnexts: %f\n", ((double)(elapsed4))/1000000000);
mexPrintf("\nj Loops: %f\n", ((double)(elapsed5))/1000000000);
mexPrintf("\nj Loops without spike prop: %f\n", ((double)(elapsed5-elapsed2))/1000000000);
t*/



/* Issue a warning if max number of spikes reached */
if(ns>=maxns)
   mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");

/* Free allocated memory */
mxFree(v);
mxFree(JnextE);
mxFree(JnextI);
mxFree(alphae);
mxFree(alphai);
mxFree(Wee);
mxFree(Wei);
mxFree(Wie);
mxFree(Wii);
mxFree(alphax);
mxFree(refstate);

}




