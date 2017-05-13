
clear
tic

rng(1)

% Number of neurons
% ignore Nx, which doesn't
% matter for this sim
Ne=10000;
Ni=10000;
Nx=1;
N=Ne+Ni;
Ne1=Ne/2;
Ni1=Ni/2;
q=Ne/N;

% Outgoing connection probabilities 
pee0=.25;
pei0=.25;
pie0=.25;
pii0=.25;

% Connection weight scalings
jee=12.5;
jei=50;
jie=20;
jii=50;

wee0=jee*pee0*q;
wei0=jei*pei0*(1-q);
wie0=jie*pie0*q;
wii0=jii*pii0*(1-q);

% Timescale, mean and variance of ffwd input
% current
taux=40;
mxe=.015*sqrt(N);
mxi=.01*sqrt(N);
vxe=0.05;
vxi=0.05;

% No ffwd input spike trains
rxe=0;
rxi=0;
jeX=0;
jiX=0;

mxe0=mxe/sqrt(N)+rxe*jeX/N;
mxi0=mxi/sqrt(N)+rxi*jiX/N;

tausyne=8;
tausyni=4;
tausynx=12;

% For balanced state to exist this vector should be decreasing
disp(sprintf('\nThis list should be decreasing for\n  a balanced state to exist: %.2f %.2f %.2f\n',mxe0/mxi0,wei0/wii0,wee0/wie0));

% and these values should be >1
disp(sprintf('\nAlso, this number should be greater than 1: %.2f\n',wii0/wee0));

disp(sprintf('\nFiring rates for large N: %.2f %.2f',1000*(mxe0*wii0-mxi0*wei0)/(wei0*wie0-wee0*wii0),1000*(mxe0*wie0-mxi0*wee0)/((wei0*wie0-wee0*wii0))))

% Neuron indices for which to record currents
nrecord=500;
Irecord=sort(randi(Ne,nrecord,1));

% Neuron params
gl=[1/15 1/10];
Cm=[1 1];
vlb=[-100 -100];
vth=[-10 -10];
DeltaT=[2 .5]; 
vT=[-50 -50]; %mV
vl=[-60 -60]; % mV
vre=[-65 -65];
tref=[1.5 .5];


% Random membrane potential initial conditions
V0min=vre(1);
V0max=vT(1);
V0=(V0max-V0min).*rand(N,1)+V0min;

% Time stuff
T=22000;
dt=.1;
Nt=round(T/dt);
Tburn=2000;
nburn=round(Tburn/dt);

% Number of outgoing connections
Kee=pee0*Ne;
Kei=pei0*Ne;
Kie=pie0*Ni;
Kii=pii0*Ni;

% Actual connection weights
Jee=(jee/sqrt(N));
Jei=-(jei/sqrt(N));
Jie=(jie/sqrt(N));
Jii=-(jii/sqrt(N));
Jex=jeX/sqrt(N);
Jix=jiX/sqrt(N);

% Build ffwd input currents
% Kx needs to be the sqrt of the desired auto-covariance
% in the sense of convolutions, ie Kx*Kx=Ax
Kx=sqrt((vxe/taux)*sqrt(2/pi))*exp(-(-6*taux:dt:6*taux).^2./((taux)^2));
sig=dt*conv(randn(Nt,1)./sqrt(dt),Kx,'same'); % white noise convolved with Kx
Ix1e=mxe+sig;
Ix1i=mxi+sig;
Ix2e=mxe+sig;
Ix2i=mxi+sig;


% Maximum number of spikes.
% Simulation will terminate with a warning if this is exceeded
maxns=N*T*.02;


[s,IF,Ie,Ii,~]=EIFTwoPopNetwork(Ix1e,Ix2e,Ix1i,Ix2i,Ne,Ni,Ne1,Ni1,Jex,Jix,Jee,Jei,Jie,Jii,rxe,rxi,Kee,Kei,Kie,Kii,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb,tausynx,tausyne,tausyni,V0,T,dt,maxns,Irecord);

s=s(:,s(1,:)>0);

reSim=1000*nnz(s(1,:)>Tburn & s(2,:)<=Ne)/(Ne*(T-Tburn))
riSim=1000*nnz(s(1,:)>Tburn & s(2,:)>Ne)/(Ni*(T-Tburn))


t0=toc;
disp(sprintf('Time for 1bcd simulation=%d sec',round(t0)))

save NetworkSimForFigure1bcd.mat

