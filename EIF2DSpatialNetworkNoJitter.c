
/* To use function from matlab, first compile by entering this into the Matlab command window:
   mex EIF2DSpatialNetwork.c
   Then call the function like this:
      [s,IF,Ie,Ii,v]=EIF2DSpatialNetworkNoJitter(sx,Nx1,Ne1,Ni1,JeF,JiF,Jee,Jei,Jie,Jii,KeF,KiF,Kee,Kei,Kie,Kii,betaeF,betaiF,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynF,tausyne,tausyni,V0,T,dt,maxns,Irecord);
 */

/* Inputs:
   sx is the feedforward presynaptic spike trains, which should be 3xNsx where Nsx is the number of spikes.
      sx(1,:) is spike times, sx(2,:) is the x-index of the neuron and sx(3,:) is the y-index
   Nx1 is the number of neurons in the feedforward layer in each direction.
      sx(2,:) and sx(3,:) should be integers between 1 and Nx1 inclusively.
   Ne1 and Ni1 are the numbers of exc and inh neurons in each direction
      So there should be Nx==Nx1^2 neurons in the ffwd layer in all, 
      Ne==Ne1^2 exc neurons in the recurrent network and Ni==Ni1^2 inh neurons
   Jab is the synaptic strength of connections from b=e,i,x to a=e,i.
   Kab is the number of projections from each cell in pop b=e,i,x to all cells in pop a=e,i
   betaab is the "width" of connections from a to b (i.e., the std of the gaussian)
      It is in units of neuron indices, so it should be proportional to Na1
   C,gl,vl,DeltaT,VT,tref,Vth,Vre,Vlb are EIF neuron params
      They are each 2x1 vectors, for exc and inh neurons separately.
      For example, C(1) is the capacitance of exc neurons and C(2) of inh neurons   
   tausynb is the time-constant of the synapses from neurons in population b=x,e,i.
     post-synaptic currents are of the form (1/tausynb)*exp(-t/tausynb) where t>0 is 
     the time evolved since the presynaptic spike.
   V0 is vector of membrane potential initial conditions.
      It should be Nx1 where N=Ne+Ni is the number of neurons in the recurrent network
      The first Ne elements are for exc neurons, the last Ni for inh neurons
   T is the total simulation time
   dt is time bin size
   maxns is maximum number of spikes allowed for all neurons together
   Irecord is a 2x(Nrecord) matrix indicating for which neurons we should record
     the synaptic inputs and membrane potential.
     Irecord(1,:) is the index in the x-direction and Irecord(2,:) in the y-direction.
     Excitatory neurons are the first Ne1 in each direction and inhibitory neurons the
     next Ni1.  For example, it Ne1=100 then Irecord(1,j)=150, Irecord(2,j)=130 means
     that the jth neuron recorded will be the inhibtiory neuron at coordinates
     (50,30)
 
   Outputs:
   s is a 3x(maxns) matrix of spikes
      s(1,:) contains spike times.
      s(2,:) contains indices of neurons that spike in the x-direction
      s(3,:) contains indices of neurons that spike in the y-direction
      Inhibitory neurons are given negative indices, exc neurons positive.
      For example the jth spike occurs at time s(1,j).
      If s(2,j)==-50 and s(3,j)==-10 then it is an inhibitory neuron spike
      at spatial coordinates (50, 10) or (50/Ni1, 10/Ni1) on the unit square.
      When there are fewer than maxns spikes, extra space in s will be filled 
      with zeros.  This should be truncated in matlab by writing s=s(:,s(1,:)>0);
   Ix,Ie,Ii,v are the recorded 
   feedforward synaptic input, exc synaptic input, inh synaptic input and voltage
   respectively.
 */



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


int printfperiod,kk,Ntref[2],Ni,Ni1,Ne2,Ni2,Kee,Kie,Kei,Kii,Kex,Kix,Ne,Ne1,j,j1,j2,k,i,N,Nt,m1,m2,maxns,ns,flagei,Nrecord,jj,nskiprecord,flag,tempflag,Nsx,Nx1,Nx;
double dt,*s,*v,*v0,*JnextE,*JnextI,*JnextX,*alphax,*alphae,*alphai,tausyne,tausyni,*alphaxr,*alphaer,*alphair,*vr,*sx,Jex,Jix,betaex,betaix,tausynx;
double Jee,Jei,Jie,Jii,betaee,betaei,betaie,betaii,T,*Irecord,*C,*Vleak,*DeltaT,*VT,*tref,*gl,*Vth,*Vre,*Vlb,xloc,yloc;
int *tempWee,*tempWei,*tempWie,*tempWii,*tempWex,*tempWix,*Wee1,*Wee2,*Wei1,*Wei2,*Wie1,*Wie2,*Wii1,*Wii2,*Wex1,*Wex2,*Wix1,*Wix2,*refstate,iXspike,jspike,postcell;
mxArray *temp1, *temp2, *temp0,*temp3,*temp4,*temp5;

/* Seed random number generator */
/* Change this if you want */
/* You could alternatively pass in a seed */
srand48(10);

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

/* Number of neurons in the ffwd layer in each direction. */ 
Nx1=(int)mxGetScalar(prhs[1]);

/* Number of exc neurons in each direction. */ 
Ne1=(int)mxGetScalar(prhs[2]);

/* Number of inh neurons in each direction. */ 
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



v0 = mxGetPr(prhs[34]);
N = mxGetM(prhs[34]);
m2 = mxGetN(prhs[34]);
if(N==1 && m2!=1)
    N=m2;

T = mxGetScalar(prhs[35]);
dt = mxGetScalar(prhs[36]);

maxns = ((int)mxGetScalar(prhs[37]));

Irecord=mxGetPr(prhs[38]);
Nrecord = mxGetN(prhs[38]);
m2 = mxGetM(prhs[38]);
if(m2!=2)
    mexErrMsgTxt("Irecord should be Nx2.");


/******
 * Finished importing variables.
 *******/

/* Total number of each type of neuron */
Ne=Ne1*Ne1;
Ni=Ni1*Ni1;
Nx=Nx1*Nx1;

/* Check for consistency with total number of neurons */
if(N!=Ne+Ni)
    mexErrMsgTxt("Ne1 and/or Ni1 not consistent with size of V0");

/* Numebr of time bins */
Nt=(int)(T/dt);


/******
 * Now allocate output variables and temporary arrays.
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

JnextX=mxMalloc(N*sizeof(double));
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

/* Inititalize variables */
for(j=0;j<N;j++){
    v[j]=v0[j]; 
    refstate[j]=0;
    JnextE[j]=0;
    JnextI[j]=0;
    alphae[j]=0;
    alphai[j]=0;
    alphax[j]=0;
    JnextX[j]=0;
}

mexPrintf("\nBuilding Network Architecure\n");
mexEvalString("drawnow;");

/* Record input currents and membrane potentials at first time bin */
for(jj=0;jj<Nrecord;jj++){
      /* Find index into local variables */
    
        /* Find index into local variables */
        j1=(int)round(Irecord[2*jj]-1);
        j2=(int)round(Irecord[2*jj+1]-1);

        /* Convert to 1D index, j */
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


      /* Store currents and membrane potentials at first time bin */
      alphaer[jj+Nrecord*0]=alphae[j];
      alphair[jj+Nrecord*0]=alphai[j];
      alphaxr[jj+Nrecord*0]=alphax[j];
      vr[jj+Nrecord*0]=v[j];
}
    


/* Initialize connections */
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
         

         }

/* Refractory states */
Ntref[0]=(int)round(tref[0]/dt);
Ntref[1]=(int)round(tref[1]/dt);

printfperiod=(int)(round(Nt/20.0));

/* Initialize number of spikes */
ns=0;

iXspike=0;

/* Print portion complete every printperiod steps */
mexPrintf("\n%d\n",printfperiod);
mexEvalString("drawnow;");

/* Main loop */
/* Exit loop and issue a warning if max number of spikes is exceeded */
for(i=1;i<Nt && ns<maxns;i++){
      
    
     /* Find all spikes in feedforward layer at this time bin */
     /* Add to corresponding elements of JnextX */     
     while(sx[iXspike*3+0]<=i*dt && iXspike<Nsx){
         
         /* Find the index of each neuron in ffwd layer that spiked */
         jspike=(int)round(sx[iXspike*3+1]-1)*Nx1+(int)round(sx[iXspike*3+2]-1);
         if(jspike<0 || jspike>=Nx){
             mexErrMsgTxt("Out of bounds index in sx.");
         }
         for(k=0;k<Kex;k++){ 
             postcell=Wex1[jspike*Kex+k]*Ne1+Wex2[jspike*Kex+k];    
             if(postcell<0 || postcell>=N)
                 mexErrMsgTxt("Out of bounds index ino JnextX");
             JnextX[postcell]+=Jex;
         }
         for(k=0;k<Kix;k++){
             postcell=Ne+Wix1[jspike*Kix+k]*Ni1+Wix2[jspike*Kix+k];
             if(postcell<0 || postcell>=N)
                 mexErrMsgTxt("Out of bounds index ino JnextX i");
             JnextX[postcell]+=Jix;
         }

         iXspike++;

     }


     for(j=0;j<N;j++){    
    
        
         /* Update synaptic variables */
         alphae[j]-=alphae[j]*(dt/tausyne);
         alphai[j]-=alphai[j]*(dt/tausyni);  
         alphax[j]-=alphax[j]*(dt/tausynx);  
      
         if(j<Ne){
             
             
             if(refstate[j]<=0)
                v[j]+=fmax((alphae[j]+alphai[j]+alphax[j]-gl[0]*(v[j]-Vleak[0])+gl[0]*DeltaT[0]*exp((v[j]-VT[0])/DeltaT[0]))*dt/C[0],Vlb[0]-v[j]);
             else{                 
                if(refstate[j]>1)
                   v[j]=Vth[0];
                else
                   v[j]=Vre[0];
                refstate[j]--;
             }
             
             
              /* If a spike occurs */
              if(v[j]>=Vth[0] && refstate[j]<=0 && ns<maxns){


                  refstate[j]=Ntref[0];
                  v[j]=Vth[0];       /* reset membrane potential */
                  s[0+3*ns]=i*dt; /* spike time */
                  s[2+3*ns]=j/Ne1+1;     /* neuron index 1 */
                  s[1+3*ns]=j%Ne1+1;     /* neuron index 2 */
                  ns++;           /* update total number of spikes */


                  /* For each postsynaptic target, propagate spike into JnextE */
                  for(k=0;k<Kee;k++)
                         JnextE[Wee1[j*Kee+k]*Ne1+Wee2[j*Kee+k]]+=Jee;
                  for(k=0;k<Kie;k++)
                         JnextE[Ne+Wie1[j*Kie+k]*Ni1+Wie2[j*Kie+k]]+=Jie;

              }
             
                 
         }
         else{ /* If cell is inhibitory */
            
             if(refstate[j]<=0)
                v[j]+=fmax((alphae[j]+alphai[j]+alphax[j]-gl[1]*(v[j]-Vleak[1])+gl[1]*DeltaT[1]*exp((v[j]-VT[1])/DeltaT[1]))*dt/C[1],Vlb[1]-v[j]);
             else{                 
                if(refstate[j]>1)
                   v[j]=Vth[1];
                else
                   v[j]=Vre[1];
                refstate[j]--;
             }
             
             
              /* If a spike occurs */
              if(v[j]>=Vth[1] && refstate[j]<=0 && ns<maxns){                                                            
                  
                  
                  refstate[j]=Ntref[1];
                  v[j]=Vth[1];       /* reset membrane potential */
                  s[0+3*ns]=i*dt; /* spike time */
                  s[2+3*ns]=-((j-Ne)/Ni1)-1;     /* neuron index 1 */
                  s[1+3*ns]=-((j-Ne)%Ni1)-1;     /* neuron index 2 */
                 
                  ns++;           /* update total number of spikes */


                  /* For each postsynaptic target, propagate spike into JnextI */
                 for(k=0;k<Kei;k++)
                        JnextI[Wei1[(j-Ne)*Kei+k]*Ne1+Wei2[(j-Ne)*Kei+k]]+=Jei;
                 for(k=0;k<Kii;k++)
                        JnextI[Ne+Wii1[(j-Ne)*Kii+k]*Ni1+Wii2[(j-Ne)*Kii+k]]+=Jii;
                  
               }
           
             
           }
          

         
        }
  

          /* Store recorded variables */
        for(jj=0;jj<Nrecord;jj++){
            
            /* Find index into local variables */
            j1=(int)round(Irecord[2*jj+0]-1);
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


              alphaer[jj+Nrecord*i]=alphae[j];
              alphair[jj+Nrecord*i]=alphai[j];
              alphaxr[jj+Nrecord*i]=alphax[j];
              vr[jj+Nrecord*i]=v[j];
              
         }    

        /* Use Jnext vectors to update synaptic variables */
        for(j=0;j<N;j++){                     
          alphae[j]+=JnextE[j]/tausyne;
          alphai[j]+=JnextI[j]/tausyni;
          alphax[j]+=JnextX[j]/tausynx;
          JnextE[j]=0;
          JnextI[j]=0; 
          JnextX[j]=0;
          
        }    
     

     /* Print percent complete every printfperiod time steps 
      * This might not actually print until the full simulation 
      * is complete due to how some versions of Matlab treat the
      * drawnow signal coming from a mex file */
     if(i%printfperiod==0){
         mexPrintf("\n%d percent complete  rate = %2.2fHz",i*100/Nt,1000*((double)(ns))/(((double)(N))*((double)(i))*dt));
         mexEvalString("drawnow;");
     }

}


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



