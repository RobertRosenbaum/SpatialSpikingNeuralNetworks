
/* To use function from matlab, first compile by entering this into the Matlab command window:
   mex LIFNetworkTimeDep0.c
   Then call the function like this:
   [s,field,alphax,alphae,alphai,v]=EIF2DFastSynFieldX(sx,Nx1,Ne1,Ni1,Jex,Jix,Jee,Jei,Jie,Jii,Kex,Kix,Kee,Kei,Kie,Kii,betaex,betaix,betaee,betaei,betaie,betaii,C,gl,vl,DeltaT,VT,tref,Vth,Vre,Vlb,tausynx,tausyne,tausyni,taujitter,V0,T,dt,maxns,Irecord,FieldLocs,FieldSigma);

 */

/* Iapp is the constant external input.  It should be a vector of size NxNxNt N is the number of cells in each direction and Nt the number of time bins in the simulation. 
   dtI is the temporal bin size used for Iapp.
   dxI is the spatial bin size used for Iapp.
   Ne is the number of excitatory cells. The first Ne indices are assumed to be excitatory and the last N-Ne are inhibitory.
   Jab is the synaptic strength of connections from b=e,i to a=e,i.  More precisely, it is the amplitude of a PSP.
   K is the number of projections from each cell.  I use a constant for simplicity, but it can be changed to a binomial random variable.
   betaab is the "width" of connections from a to b (i.e., the std of the gaussian)
   gL=1/taum, inverse of the membrane time constant
   Vth,Vre,Vlb are threshold, reset and lower boundary for LIF
   V0 is vector of membrane potential initial conditions.
   dt is bin size for time
   maxns is maximum number of spikes allowed
 
   The function returns s which is a matrix of size 2 x maxns.  
   s(1,:) contains spike times.
   s(2,:) contains indices of neurons that spike.
   When there are fewer than maxns spikes, extra space in s will be filled 
   with zeros.  This should be truncated in matlab.
 */



#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


int kk,Ntref[2],Ni,Ni1,Ne2,Ni2,Kee,Kie,Kei,Kii,Kex,Kix,Ne,Ne1,j,j1,j2,k,i,N,N1,Nt,m1,m2,maxns,ns,flagei,Nrecord,jj,nskiprecord,flag,tempflag,Nsx,Nx1,Nx,njitter,NIs,jjj;
double dt,*s,*v,*v0,*JnextE,*JnextI,*JnextX,*alphax,*alphae,*alphai,tausyne,tausyni,*alphaxr,*alphaer,*alphair,*vr,*sx,Jex,Jix,betaex,betaix,taujitter,tausynx,*Ies,*Iis,*IW;
double IMult,Jee,Jei,Jie,Jii,betaee,betaei,betaie,betaii,T,*Irecord,*C,*Vleak,*DeltaT,*VT,*tref,*gl,*Vth,*Vre,*Vlb,xloc,yloc;
int *tempWee,*tempWei,*tempWie,*tempWii,*tempWex,*tempWix,*Wee1,*Wee2,*Wei1,*Wei2,*Wie1,*Wie2,*Wii1,*Wii2,*Wex1,*Wex2,*Wix1,*Wix2,*JitterE,*JitterI,*refstate,iXspike,jspike,postcell,JitteredSpikeTime;
mxArray *temp1, *temp2, *temp0,*temp3,*temp4,*temp5;

/*t
    uint64_t        start,start1,start2,start3,start4,start5;
    uint64_t        end,end1,end2,end3,end4,end5;
    uint64_t        elapsed,elapsed1,elapsed2,elapsed3,elapsed4,elapsed5;

    elapsed=0; elapsed1=0; elapsed2=0; elapsed3=0; elapsed4=0; elapsed5=0;
t*/   

void CircRandNfun(int*,double ,double ,int ,int ,int);
/******
 * Import variables from matlab
 * This is messy looking and is specific to mex.
 * Ignore if you're implementing this outside of mex.
 *******/
sx = mxGetPr(prhs[0]);
m1 = mxGetM(prhs[0]);
Nsx = mxGetN(prhs[0]);
if(m1!=3){
    mexErrMsgTxt("sx should be Nsxx3");
}

Nx1=(int)mxGetScalar(prhs[1]);

/* Total number of neurons in each dimension. */ 
Ne1=(int)mxGetScalar(prhs[2]);

/* Total number of excitatory neurons in each dimension */
Ni1=(int)mxGetScalar(prhs[3]);

Jex=mxGetScalar(prhs[4]);
Jix=mxGetScalar(prhs[5]);

Jee= mxGetScalar(prhs[6]);
Jei= mxGetScalar(prhs[7]);
Jie= mxGetScalar(prhs[8]);
Jii= mxGetScalar(prhs[9]);

Kex=(int)mxGetScalar(prhs[10]);
Kix=(int)mxGetScalar(prhs[11]);

Kee=(int)mxGetScalar(prhs[12]);
Kei=(int)mxGetScalar(prhs[13]);
Kie=(int)mxGetScalar(prhs[14]);
Kii=(int)mxGetScalar(prhs[15]);

betaex=mxGetScalar(prhs[16]);
betaix=mxGetScalar(prhs[17]);

betaee=mxGetScalar(prhs[18]);
betaei=mxGetScalar(prhs[19]);
betaie=mxGetScalar(prhs[20]);
betaii=mxGetScalar(prhs[21]);

C=mxGetPr(prhs[22]);
m1 = mxGetN(prhs[22]);
m2 = mxGetM(prhs[22]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
gl=mxGetPr(prhs[23]);
m1 = mxGetN(prhs[23]);
m2 = mxGetM(prhs[23]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vleak=mxGetPr(prhs[24]);
m1 = mxGetN(prhs[24]);
m2 = mxGetM(prhs[24]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
DeltaT=mxGetPr(prhs[25]);
m1 = mxGetN(prhs[25]);
m2 = mxGetM(prhs[25]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
VT=mxGetPr(prhs[26]);
m1 = mxGetN(prhs[26]);
m2 = mxGetM(prhs[26]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
tref=mxGetPr(prhs[27]);
m1 = mxGetN(prhs[27]);
m2 = mxGetM(prhs[27]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vth=mxGetPr(prhs[28]);
m1 = mxGetN(prhs[28]);
m2 = mxGetM(prhs[28]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vre=mxGetPr(prhs[29]);
m1 = mxGetN(prhs[29]);
m2 = mxGetM(prhs[29]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vlb=mxGetPr(prhs[30]);
m1 = mxGetN(prhs[30]);
m2 = mxGetM(prhs[30]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");


tausynx=mxGetScalar(prhs[31]);
tausyne=mxGetScalar(prhs[32]);
tausyni=mxGetScalar(prhs[33]);

taujitter=mxGetScalar(prhs[34]);


v0 = mxGetPr(prhs[35]);
N = mxGetM(prhs[35]);
m2 = mxGetN(prhs[35]);
if(N==1 && m2!=1)
    N=m2;

T = mxGetScalar(prhs[36]);
dt = mxGetScalar(prhs[37]);

maxns = ((int)mxGetScalar(prhs[38]));

Irecord=mxGetPr(prhs[39]);
Nrecord = mxGetN(prhs[39]);
m2 = mxGetM(prhs[39]);
if(m2!=2)
    mexErrMsgTxt("Irecord should be Nx2.");


Nt=(int)(T/dt);

Ies=mxGetPr(prhs[40]);
NIs = mxGetM(prhs[40]);
m2 = mxGetN(prhs[40]);
if(m2!=Nt)
    mexErrMsgTxt("Ie should be ....");
mexPrintf("\nNIs=%d\n",NIs);


Iis=mxGetPr(prhs[41]);
m1 = mxGetM(prhs[41]);
m2 = mxGetN(prhs[41]);
if(m2!=Nt || m1!=NIs)
    mexErrMsgTxt("Ii should be same as Ie.");

IW=mxGetPr(prhs[42]);
m1 = mxGetM(prhs[42]);
m2 = mxGetN(prhs[42]);
if(m1!=NIs||m2!=N)
    mexErrMsgTxt("IW should be ...");
    
mexPrintf("\n%d %d\n",NIs,Nt);
mexEvalString("drawnow;");

/******
 * Finished importing variables.
 *******/

for (int krand = 0; krand < 40; krand++) {
	mexPrintf("%f\n", drand48());
}

Ne=Ne1*Ne1;
Ni=Ni1*Ni1;
N1=Ne1+Ni1;
Nx=Nx1*Nx1;

if(N!=Ne+Ni)
    mexErrMsgTxt("Ne1 and/or Ni1 not consistent with size of V0");


/* Number of bins for jitter storage */
njitter=(int)round(taujitter*5/dt);

/******
 * Now allocate new variables.
 * This is also mex specific.  Use malloc in C, etc.
 *****/

/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(3, maxns, mxREAL);
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
temp0=mxCreateDoubleMatrix(N, 1, mxREAL);
v = mxGetPr(temp0);

temp1=mxCreateDoubleMatrix(N, 1, mxREAL);
JnextE = mxGetPr(temp1);

temp2=mxCreateDoubleMatrix(N, 1, mxREAL);
JnextI = mxGetPr(temp2);

temp3=mxCreateDoubleMatrix(N, 1, mxREAL);
alphae = mxGetPr(temp3);

temp4=mxCreateDoubleMatrix(N, 1, mxREAL);
alphai = mxGetPr(temp4);

JnextX=mxMalloc(N*njitter*sizeof(double));
alphax=mxMalloc(N*sizeof(double));

/* Vectors for postsynaptic connections */
Wee1=mxMalloc(Ne*Kee*sizeof(int));
Wee2=mxMalloc(Ne*Kee*sizeof(int));
Wei1=mxMalloc(Ni*Kei*sizeof(int));
Wei2=mxMalloc(Ni*Kei*sizeof(int));
Wie1=mxMalloc(Ne*Kie*sizeof(int));
Wie2=mxMalloc(Ne*Kie*sizeof(int));
Wii1=mxMalloc(Ni*Kii*sizeof(int));
Wii2=mxMalloc(Ni*Kii*sizeof(int));
Wex1=mxMalloc(Nx*Kex*sizeof(int));
Wex2=mxMalloc(Nx*Kex*sizeof(int));
Wix1=mxMalloc(Nx*Kix*sizeof(int));
Wix2=mxMalloc(Nx*Kix*sizeof(int));
/*JitterE=mxMalloc(Nx*Kex*sizeof(int));
JitterI=mxMalloc(Nx*Kix*sizeof(int));*/
tempWee=mxMalloc(Kee*sizeof(int));
tempWei=mxMalloc(Kei*sizeof(int));
tempWie=mxMalloc(Kie*sizeof(int));
tempWii=mxMalloc(Kii*sizeof(int));
tempWex=mxMalloc(Kex*sizeof(int));
tempWix=mxMalloc(Kix*sizeof(int));
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
    for(i=0;i<njitter;i++)
        JnextX[j*njitter+i]=0;
}


for(jj=0;jj<Nrecord;jj++){
      /* Find index into local variables */
    
        /* Find index into local variables */
        j1=(int)round(Irecord[2*jj]-1);
        j2=(int)round(Irecord[2*jj+1]-1);

        if(j1<Ne1 && j2<Ne1){                 
            j=j1+Ne1*j2;
        }
        else
          if(j1>=Ne1 && j2>=Ne1){
            j=(j1-Ne1)+(j2-Ne1)*Ni1+Ne;
          }
          else
            mexErrMsgTxt("Indices in Irecord must have both terms <Ne1 or both terms >Ne1");


                
      if(j>=N || j<0){         
         mexErrMsgTxt("Irecord contains out of bounds indices.");
      }

      alphaer[jj+Nrecord*0]=alphae[j];
      alphair[jj+Nrecord*0]=alphai[j];
      alphaxr[jj+Nrecord*0]=alphax[j];
      vr[jj+Nrecord*0]=v[j];
}
    

srand48(10);


/*t
start = mach_absolute_time();
t*/
/* Initialize connection matrix */
for(j=0;j<Ne;j++){
                
         /* Find index of cell along each dimension */
         j1=j/Ne1;
         j2=j%Ne1;
         
         /* Generate vectors of exc and inh postsynaptic targets */
         CircRandNfun(tempWee,(double)j1,betaee,0,Ne1-1,Kee);
         for(k=0;k<Kee;k++)
             Wee1[j*Kee+k]=tempWee[k];
         CircRandNfun(tempWee,(double)j2,betaee,0,Ne1-1,Kee);
         for(k=0;k<Kee;k++)
             Wee2[j*Kee+k]=tempWee[k];      
         
         CircRandNfun(tempWie,((double)j1*((double)Ni1/(double)Ne1)),betaie,0,Ni1-1,Kie);
         for(k=0;k<Kie;k++)
             Wie1[j*Kie+k]=tempWie[k];         
         CircRandNfun(tempWie,((double)j2*((double)Ni1/(double)Ne1)),betaie,0,Ni1-1,Kie);
         for(k=0;k<Kie;k++)
             Wie2[j*Kie+k]=tempWie[k];                   
     }


for(j=0;j<Ni;j++){
                         
         /* Find index of cell along each dimension */
         j1=j/Ni1;
         j2=j%Ni1;
         
        /* Generate vectors of exc and inh postsynaptic targets */        
        CircRandNfun(tempWei,((double)j1*((double)Ne1/(double)Ni1)),betaei,0,Ne1-1,Kei);
        for(k=0;k<Kei;k++)
            Wei1[j*Kei+k]=tempWei[k];        
        CircRandNfun(tempWei,((double)j2*((double)Ne1/(double)Ni1)),betaei,0,Ne1-1,Kei);
        for(k=0;k<Kei;k++)
            Wei2[j*Kei+k]=tempWei[k];      
        
        CircRandNfun(tempWii,(double)j1,betaii,0,Ni1-1,Kii);
        for(k=0;k<Kii;k++)
            Wii1[j*Kii+k]=tempWii[k];           
        CircRandNfun(tempWii,(double)j2,betaii,0,Ni1-1,Kii);
        for(k=0;k<Kii;k++)
            Wii2[j*Kii+k]=tempWii[k];    
     }


/* Initialize connection matrix */
for(j=0;j<Nx;j++){
                
         /* Find index of cell along each dimension */
         j1=j/Nx1;
         j2=j%Nx1;
         
         /* Generate vectors of exc and inh postsynaptic targets */
         CircRandNfun(tempWex,((double)j1*((double)Ne1/(double)Nx1)),betaex,0,Ne1-1,Kex);
         for(k=0;k<Kex;k++)
             Wex1[j*Kex+k]=tempWex[k];
         CircRandNfun(tempWex,((double)j2*((double)Ne1/(double)Nx1)),betaex,0,Ne1-1,Kex);
         for(k=0;k<Kex;k++)
             Wex2[j*Kex+k]=tempWex[k];      
         
         CircRandNfun(tempWix,((double)j1*((double)Ni1/(double)Nx1)),betaix,0,Ni1-1,Kix);
         for(k=0;k<Kix;k++)
             Wix1[j*Kix+k]=tempWix[k];
         CircRandNfun(tempWix,((double)j2*((double)Ni1/(double)Nx1)),betaix,0,Ni1-1,Kix);
         for(k=0;k<Kix;k++)
             Wix2[j*Kix+k]=tempWix[k];    
         
         
         
/*         jspike=(int)round(sx[iXspike*3+1]-1)*Nx1+(int)round(sx[iXspike*3+2]-1);*/
/*         if(jspike<0 || jspike>=Nx)
             mexErrMsgTxt("Out of bounds index in sx.");*/
/*         CircRandNfun(tempWex,(((double)njitter)/2),taujitter/dt,0,njitter-1,Kex,(double)(j+4*Ne+4*Ni+4*Nx));  
         for(k=0;k<Kex;k++)
             JitterE[j*Kex+k]=tempWex[j];
         CircRandNfun(tempWix,(((double)njitter)/2),taujitter/dt,0,njitter-1,Kix,(double)(j+4*Ne+4*Ni+5*Nx));  
         for(k=0;k<Kix;k++)
             JitterI[j*Kix+k]=tempWix[j];
*/         
         
         }

/*t
end = mach_absolute_time();
elapsed = end - start;
mexPrintf("\nBuilding Ws: %" PRIu64 "\n", elapsed);
t*/


/*
mexPrintf("\n");
for(j=0;j<Ne;j++)
for(k=0;k<Kie;k++)
    if(Wie1[j*Kie+k]==100 && Wie2[j*Kie+k]==100)
        mexPrintf("%d\n",j);

mexPrintf("\n\n Break \n\n");

mexPrintf("\n");
for(j=0;j<Ne;j++)
for(k=0;k<Kie;k++)
    if(Wie1[j*Kie+k]==100 && Wie2[j*Kie+k]==101)
        mexPrintf("%d\n",j);
*/

Ntref[0]=(int)round(tref[0]/dt);
Ntref[1]=(int)round(tref[1]/dt);


/*mexPrintf("\n%d %f\n",Nt,Iappe[100]);*/

/* Initialize number of spikes */
ns=0;

flag=0;

iXspike=0;

srand48((double)(4*Ne+4*Ni+4*Nx));

mexPrintf("\nnjitter=%d,taujitter=%f,dt=%f\n",njitter,taujitter,dt);

/*t start = mach_absolute_time(); t*/
/* Time loop */
/* Exit loop and issue a warning if max number of spikes is exceeded */
for(i=1;i<Nt && ns<maxns;i++){
      
      while(sx[iXspike*3+0]<=i*dt && iXspike<Nsx){
         
         
         jspike=(int)round(sx[iXspike*3+1]-1)*Nx1+(int)round(sx[iXspike*3+2]-1);
         if(jspike<0 || jspike>=Nx){
             mexPrintf("\n %d %d %d %d %d %d\n",(int)round(sx[iXspike*3+0]/dt),iXspike,i,jspike,(int)round(sx[iXspike*3+1]-1),(int)round(sx[iXspike*3+2]-1));
             mexErrMsgTxt("Out of bounds index in sx.");
         }
         CircRandNfun(tempWex,(double)((njitter/2)),taujitter/dt,0,njitter-1,Kex);
         for(k=0;k<Kex;k++){ 
             if(jspike*Kex+k>=Nx*Kex || jspike*Kex+k<0)
                 mexErrMsgTxt("Out of bounds jspike E.");
             postcell=Wex1[jspike*Kex+k]*Ne1+Wex2[jspike*Kex+k];    
             JitteredSpikeTime=tempWex[k]+i;
             if(postcell*njitter+JitteredSpikeTime%njitter<0 || postcell*njitter+JitteredSpikeTime%njitter>=N*njitter)
                 mexErrMsgTxt("Out of bounds index ino JnextX");
             JnextX[postcell*njitter+JitteredSpikeTime%njitter]+=Jex;
         }
         CircRandNfun(tempWix,(double)((njitter/2)),taujitter/dt,0,njitter-1,Kix);
         for(k=0;k<Kix;k++){
             if(jspike*Kix+k>=Nx*Kix || jspike*Kix+k<0)
                 mexErrMsgTxt("Out of bounds jspike I.");
             postcell=Ne+Wix1[jspike*Kix+k]*Ni1+Wix2[jspike*Kix+k];
             JitteredSpikeTime=tempWix[k]+i;
             if(postcell*njitter+JitteredSpikeTime%njitter<0 || postcell*njitter+JitteredSpikeTime%njitter>=N*njitter)
                 mexErrMsgTxt("Out of bounds index ino JnextX i");
             JnextX[postcell*njitter+JitteredSpikeTime%njitter]+=Jix;
         }

         iXspike++;

     }

	  for (j = 0; j<N; j++){


		  /* Update synaptic variables */
		  alphae[j] -= alphae[j] * (dt / tausyne);
		  alphai[j] -= alphai[j] * (dt / tausyni);
		  alphax[j] -= alphax[j] * (dt / tausynx);

		  if (j<Ne){


			  if (refstate[j] <= 0)
				  v[j] += fmax((alphae[j] + alphai[j] + (1+Ies[i]*IW[j])*alphax[j] - gl[0] * (v[j] - Vleak[0]) + gl[0] * DeltaT[0] * exp((v[j] - VT[0]) / DeltaT[0]))*dt / C[0], Vlb[0] - v[j]);
			  else{
				  if (refstate[j]>1)
					  v[j] = Vth[0];
				  else
					  v[j] = Vre[0];
				  refstate[j]--;
			  }


			  /* If a spike occurs */
			  if (v[j] >= Vth[0] && refstate[j] <= 0 && ns<maxns){


				  refstate[j] = Ntref[0];
				  v[j] = Vth[0];       /* reset membrane potential */
				  s[0 + 3 * ns] = i*dt; /* spike time */
				  s[2 + 3 * ns] = j / Ne1 + 1;     /* neuron index 1 */
				  s[1 + 3 * ns] = j%Ne1 + 1;     /* neuron index 2 */
				  ns++;           /* update total number of spikes */


				  /* For each postsynaptic target, propagate spike into JnextE */
				  for (k = 0; k<Kee; k++)
					  JnextE[Wee1[j*Kee + k] * Ne1 + Wee2[j*Kee + k]] += Jee;
				  for (k = 0; k<Kie; k++)
					  JnextE[Ne + Wie1[j*Kie + k] * Ni1 + Wie2[j*Kie + k]] += Jie;

			  }


		  }
		  else{ /* If cell is inhibitory */

			  if (refstate[j] <= 0)
				  v[j] += fmax((alphae[j] + alphai[j] + (1+Iis[i]*IW[j])*alphax[j] - gl[1] * (v[j] - Vleak[1]) + gl[1] * DeltaT[1] * exp((v[j] - VT[1]) / DeltaT[1]))*dt / C[1], Vlb[1] - v[j]);
			  else{
				  if (refstate[j]>1)
					  v[j] = Vth[1];
				  else
					  v[j] = Vre[1];
				  refstate[j]--;
			  }


			  /* If a spike occurs */
			  if (v[j] >= Vth[1] && refstate[j] <= 0 && ns<maxns){


				  refstate[j] = Ntref[1];
				  v[j] = Vth[1];       /* reset membrane potential */
				  s[0 + 3 * ns] = i*dt; /* spike time */
				  s[2 + 3 * ns] = -((j - Ne) / Ni1) - 1;     /* neuron index 1 */
				  s[1 + 3 * ns] = -((j - Ne) % Ni1) - 1;     /* neuron index 2 */

				  ns++;           /* update total number of spikes */


				  /* For each postsynaptic target, propagate spike into JnextI */
				  for (k = 0; k<Kei; k++)
					  JnextI[Wei1[(j - Ne)*Kei + k] * Ne1 + Wei2[(j - Ne)*Kei + k]] += Jei;
				  for (k = 0; k<Kii; k++)
					  JnextI[Ne + Wii1[(j - Ne)*Kii + k] * Ni1 + Wii2[(j - Ne)*Kii + k]] += Jii;

			  }


		  }



	  }


	  /* Store recorded variables */
	  for (jj = 0; jj<Nrecord; jj++){

		  /* Find index into local variables */
		  j1 = (int)round(Irecord[2 * jj + 0] - 1);
		  j2 = (int)round(Irecord[2 * jj + 1] - 1);

		  if (j1<Ne1 && j2<Ne1){
			  j = j1 + Ne1*j2;
		  }
		  else
			  if (j1 >= Ne1 && j2 >= Ne1){
				  j = (j1 - Ne1) + (j2 - Ne1)*Ni1 + Ne;
			  }
			  else
				  mexErrMsgTxt("Indices in Irecord must have both terms <Ne1 or both terms >Ne1");


		  alphaer[jj + Nrecord*i] = alphae[j];
		  alphair[jj + Nrecord*i] = alphai[j];
		  alphaxr[jj + Nrecord*i] = alphax[j];
		  vr[jj + Nrecord*i] = v[j];

	  }

	  /* Use Jnext vectors to update synaptic variables */
	  for (j = 0; j<N; j++){
		  alphae[j] += JnextE[j] / tausyne;
		  alphai[j] += JnextI[j] / tausyni;
		  alphax[j] += JnextX[j*njitter + (i%njitter)] / tausynx;
		  JnextE[j] = 0;
		  JnextI[j] = 0;
		  JnextX[j*njitter + (i%njitter)] = 0;

	  }



  
    
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
mxDestroyArray(temp0);
mxDestroyArray(temp1);
mxDestroyArray(temp2);
mxDestroyArray(temp3);
mxDestroyArray(temp4);
mxFree(Wee1);
mxFree(Wee2);
mxFree(Wei1);
mxFree(Wei2);
mxFree(Wie1);
mxFree(Wie2);
mxFree(Wii1);
mxFree(Wii2);
mxFree(Wex1);
mxFree(Wex2);
mxFree(Wix1);
mxFree(Wix2);
mxFree(tempWee);
mxFree(tempWei);
mxFree(tempWie);
mxFree(tempWii);
mxFree(tempWix);
mxFree(tempWex);
mxFree(JnextX);
mxFree(alphax);
mxFree(refstate);

}




/* This function generates n random variables that are integers between min and max.
 * The distribution of these random variables is like a Gaussian distribution with
 * mean mu and std sigma, but rounded to the nearest integer and wrapped around the 
 * interval [min,max]. The values are stored in the location pointed to by z, so you 
 * better make sure that you allocated room for at least n integers in z.
 */
void CircRandNfun(int* z, double mu,double sigma,int min,int max,int n){
int i;
double u1,u2;
int matlabmod(int, int);

    for(i=1;i<n;i+=2){    
        u1=drand48();
        u2=drand48();

         z[i-1]=matlabmod((int)round(sigma*sqrt(-2*log(u1))*cos(2*M_PI*u2)+mu)-min,max-min+1)+min;
         z[i]=matlabmod((int)round(sigma*sqrt(-2*log(u1))*sin(2*M_PI*u2)+mu)-min,max-min+1)+min;
    }

    if(i==n){
        u1=drand48();
        u2=drand48();
        z[i-1]=(int)round(sigma*sqrt(-2*log(u1))*cos(2*M_PI*u2)+mu);
        z[i-1]=matlabmod(z[i-1]-min,max-min+1)+min;
    }
    
}


/* Implements a "mod" function that behaves like matlab's version
 * of mod instead of like % in C.  I forgot what the difference is, 
 * but remember that it's important, especially when a is negative.
 */
int matlabmod(int a, int b){
int c;
    c=a%b;
while(c<0)
    c+=b;
return c;
}



