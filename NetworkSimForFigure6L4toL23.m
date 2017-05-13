
clear

load NetworkSimForFigure6LGNtoL4;
sx=s0(:,s0(2,:)>0);
clear s0;

rng(1)

tic

% Number of neurons in network
Ne1=200; % Number in each direction
Ni1=100;
Nx1=max(max(sx(2,:)),max(sx(3,:)));
Ne=Ne1*Ne1; % Total numbers
Ni=Ni1*Ni1;
Nx=Nx1*Nx1;
N1=Ne1+Ni1;
N=Ne+Ni;

q=Ne/N;  % Exc-Inh ratio

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

jeX=220;
jiX=220;
rX=reSim;
pex0=.25*(75^2/Nx);
pix0=.08*(75^2/Nx);

wex0=jeX*pex0*(Nx/N);
wix0=jiX*pix0*(Nx/N);

sigmaeX=.05;
sigmaiX=.05;
sigmaRe=.15;
sigmaRi=.03;
sigmac=.25;


numfactors=1;

taujitter=2;
tauI=40;
strengthI=.5;
tausyne=6;
tausyni=5;
tausynx=tausyne;

nrecord0=2%200;
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


% Build timescale and weights
% for gain modulation
tauKx=-5*tauI:dt:5*tauI;
Kx=exp(-tauKx.^2./(2*tauI^2));
Kx=Kx./sum(Kx);

Ies=zeros(numfactors,Nt);
Iis=zeros(numfactors,Nt);

for j=1:numfactors
    I=randn(Nt,1);
    I=conv(I,Kx,'same');
    I=I./std(I);
    Ies(j,:)=I*strengthI;
    Iis(j,:)=I*strengthI;
end

[x,y]=meshgrid((1:Ne1)/Ne1,(1:Ne1)/Ne1);
IWe=WrappedGauss(x(:),0.5,sigmac,5).*WrappedGauss(y(:),0.5,sigmac,5);
[x,y]=meshgrid((1:Ni1)/Ni1,(1:Ni1)/Ni1);
IWi=WrappedGauss(x(:),0.5,sigmac,5).*WrappedGauss(y(:),0.5,sigmac,5);
IW(1,:)=[IWe' IWi'];



% Maximum number of spikes.
% Simulation will terminate with a warning if this is exceeded
maxns=N*T*.1;

% Random initial membrane potentials
V0=(V0max-V0min).*rand(N,1)+V0min;

[s,~,~,~,~]=EIF2DSpatialNetworkGainMod(sx,Nx1,Ne1,Ni1,Jex,Jix,Jee,Jei,Jie,Jii,Kex,Kix,Kee,Kei,Kie,Kii,betaex,betaix,betaee,betaei,betaie,betaii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynx,tausyne,tausyni,taujitter,V0,T,dt,maxns,Irecord,Ies,Iis,IW);    
s=s(:,s(2,:)~=0);

reSim=nnz(s(1,:)>Tburn & s(2,:)>0)/(Ne*(T-Tburn))
riSim=nnz(s(1,:)>Tburn & s(2,:)<0)/(Ni*(T-Tburn))


t0=toc

disp(sprintf('Time for Fig6L4toL23 simulation=%d sec',round(t0)))

save NetworkSimForFigure6L4toL23;

