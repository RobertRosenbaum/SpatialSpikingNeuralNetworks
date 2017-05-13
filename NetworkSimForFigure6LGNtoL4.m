
clear

rng(1);

tic

% The comments and variable names in this
% code are not as neat and readable as the
% other sims.  See the sims for Figure3 and
% Figure4 for cleaner versions of similar sims

% Number of neurons in network
Ne1=200;
Ni1=100;
Nx1=75;
Ne=Ne1*Ne1;
Ni=Ni1*Ni1;
Nx=Nx1*Nx1;
N1=Ne1+Ni1;
N=Ne+Ni;

q=Ne/N;

% Scaling factors for connection probability (k) and connection strength
% (j)
pee0=.05;
pei0=.05;

jee=40;
jei=400;

% To i1
pie0=.05;
pii0=.05;

jie=120; 
jii=400;

wee0=jee*pee0*q;
wei0=jei*pei0*(1-q);
wie0=jie*pie0*q;
wii0=jii*pii0*(1-q);

jeX=120;
jiX=120;
rX=.005;
pex0=.25;
pix0=.08;

wex0=jeX*pex0*(Nx/N);
wix0=jiX*pix0*(Nx/N);

sigmaeX=.1;
sigmaiX=.1;
sigmaRe=.05;
sigmaRi=.03;

taujitter=2;
tausyne=6;
tausyni=5;
tausynx=tausyne;

nrecord0=2%00;
Irecord=randi(Ne1,2,nrecord0);


% For balanced state to exist this vector should be decreasing
disp(sprintf('\nThis list should be decreasing for\n  a balanced state to exist: %.2f %.2f %.2f\n',wex0/wix0,wei0/wii0,wee0/wie0));

% and these values should be >1
disp(sprintf('\nAlso, this number should be greater than 1: %.2f\n',wii0/wee0));

disp(sprintf('\nFiring rates for large N: %.2f %.2f',rX*1000*(wex0*wii0-wix0*wei0)/(wei0*wie0-wee0*wii0),rX*1000*(wex0*wie0-wix0*wee0)/(wei0*wie0-wee0*wii0)))


gl=[1/15 1/10];
Cm=[1 1];
vlb=[-100 -100];
vth=[-10 -10];
DeltaT=[2 .5];
vT=[-50 -50]; %mV
vl=[-60 -60]; % mV
vre=[-65 -65];
tref=[1.5 .5];


V0min=vre(1);
V0max=vT(1);


% Width of recurrent connections
sigmaee=sigmaRe;
sigmaie=sigmaRe;
sigmaei=sigmaRi;
sigmaii=sigmaRi;


%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%

T=22000;
dt=.1;
Nt=round(T/dt);
Tburn=2000;
nburn=round(Tburn/dt);

betaee=sigmaee*(Ne1);
betaei=sigmaei*(Ne1);
betaie=sigmaie*(Ni1);
betaii=sigmaii*(Ni1);
betaex=sigmaeX*(Ne1);
betaix=sigmaiX*(Ni1);

Kee=pee0*Ne;
Kei=pei0*Ne;

Kie=pie0*Ni;
Kii=pii0*Ni;

Kex=pex0*Ne;
Kix=pix0*Ni;

Jee=(jee/sqrt(N));
Jei=-(jei/sqrt(N));

Jie=(jie/sqrt(N));
Jii=-(jii/sqrt(N));

Jex=jeX/sqrt(N);
Jix=jiX/sqrt(N);

% Generate spike times for the ffwd network
tempspikes=sort(T*rand(poissrnd(rX*Nx*T),1));

% Store in the data format expected by the C code,
% where neuron indices are assigned randomly
sx=zeros(3,numel(tempspikes));
sx(1,:)=tempspikes;
clear tempspikes;
sx(2,:)=ceil(rand(1,size(sx,2))*Nx1);
sx(3,:)=ceil(rand(1,size(sx,2))*Nx1);


% Maximum number of spikes.
% Simulation will terminate with a warning if this is exceeded
maxns=N*T*.1;

% Random initial membrane potentials
V0=(V0max-V0min).*rand(N,1)+V0min;


% Zero gain modulation
% Only use non-zero modulation in L23
% No modulation in L4 (here)
Ies=zeros(1,Nt);
Iis=zeros(1,Nt);
[x,y]=meshgrid((1:Ne1)/Ne1,(1:Ne1)/Ne1);
IWe=zeros(Ne,1);
IWi=zeros(Ni,1);
IW(1,:)=[IWe' IWi'];


% Simulate Network
[s0,~,~,~,~]=EIF2DSpatialNetworkGainMod(sx,Nx1,Ne1,Ni1,Jex,Jix,Jee,Jei,Jie,Jii,Kex,Kix,Kee,Kei,Kie,Kii,betaex,betaix,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynx,tausyne,tausyni,taujitter,V0,T,dt,maxns,Irecord,Ies,Iis,IW);    
s0=s0(:,s0(2,:)~=0);

reSim=nnz(s0(1,:)>Tburn & s0(2,:)>0)/(Ne*(T-Tburn))
riSim=nnz(s0(1,:)>Tburn & s0(2,:)<0)/(Ni*(T-Tburn))

t0=toc

disp(sprintf('Time for Fig6LGNtoL4 simulation=%d sec',round(t0)))

save NetworkSimForFigure6LGNtoL4.mat

