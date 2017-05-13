
clear
close all

figure

rng(1)

%%%%%%%%%%%%%
% First panels b,c,d
%%%%%%%%%%%%%
load NetworkSimForFigure1bcd.mat;


%%%%%%%%%%%%%%%%%%%%%
% Make raster plot
%%%%%%%%%%%%%%%%%%%%%

% Size of circles plotted
mrkrsz=5;

% Time window over which to plot in ms
Tmax=2000;
Tmin=1000;

% Number of cells to plot
nplot=500;

% indices to plot
Iplot=find(s(1,:)>=Tmin & s(1,:)<=Tmax & s(2,:)<=nplot);

subplot(2,4,2)
plot((s(1,Iplot))/1000,s(2,Iplot),'k.','MarkerSize',mrkrsz)
axis tight

%%%%%%%%%%%%%%%%%%%%%
% Plot spike count correlation distribution
%%%%%%%%%%%%%%%%%%%%%

% counting window size
winsize=250;

% number of cells to sample
nc=1000;

% Edges for histogram
edgest=0:winsize:T;
edgesi=(1:Ne+1)-.01;
edges={edgest,edgesi};

% Find excitatory spikes,
% store into s0, which has the structure
% needed for hist3
Is=find(s(2,:)>0);
s0=s(:,Is)';

% Get 2D histogram of spike indices and times
counts=hist3(s0,'Edges',edges);

% Get rid of edges, 
% the last element in each
% direction is zero
counts=counts(ceil(Tburn/winsize):end-1,1:end-1);

% Find neurons with rates >=2Hz
Igood=find(mean(counts)/winsize>1/1000);

% Randomly choose nc cells to use
% their indices are stored in Inds
temp=randperm(numel(Igood),nc);
Inds=Igood(temp);

% Store their spike counts
counts=counts(:,Inds);

% Compute their pairwise correlations
C=corrcoef(counts);

% Only keep the lower-half of the correlation matrix
% store into SCorree
[II,JJ]=meshgrid(1:size(C,1),1:size(C,1));
SCorree=C(II>JJ & isfinite(C));

% Plot correlation distribution
[h,b]=hist(SCorree,100);
h=h./trapz(b,h);
subplot(2,4,5)
plot(b,h,'k','Linewidth',2)
temp=axis;
temp(1)=-1;
temp(2)=1;
axis(temp);


%%%%%%%%%%%%%
% Plot currents
%%%%%%%%%%%%%

% Time window over which to plot in ms
Tmax=2000;
Tmin=1000;

% Time constant of filter
tauKc=15;

% Low-pass filter
Kc=exp(-abs(-5*tauKc:dt:5*tauKc)/tauKc);
Kc=Kc/sum(Kc);

IF0=Ix1e';%mean(IF);
Ie0=mean(Ie);
Ii0=mean(Ii);

% Low-pass filter currents
IF0=conv2(IF0,Kc,'same');
Ie0=conv2(Ie0,Kc,'same');
Ii0=conv2(Ii0,Kc,'same');

% Compute neurons' rehobases
Grheobase=CalcRheoBaseEIF(Cm(1),gl(1),vl(1),DeltaT(1),vT(1),vth(1),vl(1),0,10,1000,dt,20);

subplot(2,4,6)
time=dt:dt:Nt;
I=find(time>=Tmin & time<=Tmax);
plot(time(I),(Ie0(I)+Ii0(I)-mean(Ie0(I)+Ii0(I)))/Grheobase,'r','LineWidth',1.5)
hold on
plot(time(I),(IF0(I)-mean(IF0(I)))/Grheobase,'b','LineWidth',1.5)
plot(time(I),(Ie0(I)+Ii0(I)+IF0(I)-mean(Ie0(I)+Ii0(I)+IF0(I)))/Grheobase,'k','LineWidth',1.5)


%%%%%%%%%%%%%
% Now panels f,g,h
%%%%%%%%%%%%%
load NetworkSimForFigure1fgh.mat;


%%%%%%%%%%%%%%%%%%%%%
% Make raster plot
%%%%%%%%%%%%%%%%%%%%%

% Size of circles plotted
mrkrsz=5;

% Time window over which to plot in ms
Tmax=2000;
Tmin=1000;

% Number of cells to plot
nplot=500;

% indices to plot
Iplot1=find(s(1,:)>=Tmin & s(1,:)<=Tmax & s(2,:)<=nplot/2);
Iplot2=find(s(1,:)>=Tmin & s(1,:)<=Tmax & s(2,:)>Ne/2 & s(2,:)<=Ne/2+nplot/2);
Iplot=[Iplot1 Iplot2];
neuroninds=s(2,Iplot);
neuroninds(neuroninds>Ne/2)=neuroninds(neuroninds>Ne/2)-Ne/2+nplot/2;

subplot(2,4,4)
plot((s(1,Iplot))/1000,neuroninds,'k.','MarkerSize',mrkrsz)
axis tight

%%%%%%%%%%%%%%%%%%%%%
% Plot spike count correlation distribution
%%%%%%%%%%%%%%%%%%%%%

% counting window size
winsize=250;

% number of cells to sample
nc=1000;

% Edges for histogram
edgest=0:winsize:T;
edgesi=(1:Ne+1)-.01;
edges={edgest,edgesi};

% Find excitatory spikes,
% store into s0, which has the structure
% needed for hist3
Is=find(s(2,:)>0);
s0=s(:,Is)';

% Get 2D histogram of spike indices and times
counts=hist3(s0,'Edges',edges);

% Get rid of edges, 
% the last element in each
% direction is zero
counts=counts(ceil(Tburn/winsize):end-1,1:end-1);

% Find neurons with rates >=2Hz
Igood=find(mean(counts)/winsize>1/1000);

% Randomly choose nc cells to use
% their indices are stored in Inds
temp=randperm(numel(Igood),nc);
Inds=Igood(temp);

% Store their spike counts
counts=counts(:,Inds);

% Compute their pairwise correlations
C=corrcoef(counts);

% Compute their "distances"
pops=(Inds<=(Ne/2));
[pops1,pops2]=meshgrid(pops,pops);
d=abs(pops1-pops2);


% Only keep the lower-half of the correlation matrix
% store into SCorree
[II,JJ]=meshgrid(1:size(C,1),1:size(C,1));
SCorree=C(II>JJ & isfinite(C));
d=d(II>JJ & isfinite(C));

% Plot correlation distribution
subplot(2,4,7)
[h,b]=hist(SCorree,100);
h=h./trapz(b,h);
plot(b,h,'k','Linewidth',2)
hold on
[h,b]=hist(SCorree(d==0),100);
h=h./trapz(b,h);
plot(b,h,'m','Linewidth',2)
[h,b]=hist(SCorree(d==1),100);
h=h./trapz(b,h);
plot(b,h,'g','Linewidth',2)
temp=axis;
temp(1)=-1;
temp(2)=1;
axis(temp);


%%%%%%%%%%%%%
% Plot currents
%%%%%%%%%%%%%

% Time window over which to plot in ms
Tmax=2000;
Tmin=1000;

% Time constant of filter
tauKc=15;

% Low-pass filter
Kc=exp(-abs(-5*tauKc:dt:5*tauKc)/tauKc);
Kc=Kc/sum(Kc);

IF0=Ix1e';%mean(IF);
Ie0=mean(Ie);
Ii0=mean(Ii);

% Low-pass filter currents
IF0=conv2(IF0,Kc,'same');
Ie0=conv2(Ie0,Kc,'same');
Ii0=conv2(Ii0,Kc,'same');

% Compute neurons' rehobases
Grheobase=CalcRheoBaseEIF(Cm(1),gl(1),vl(1),DeltaT(1),vT(1),vth(1),vl(1),0,10,1000,dt,20);

subplot(2,4,8)
time=dt:dt:Nt;
I=find(time>=Tmin & time<=Tmax);
plot(time(I),(Ie0(I)+Ii0(I)-mean(Ie0(I)+Ii0(I)))/Grheobase,'r','LineWidth',1.5)
hold on
plot(time(I),(IF0(I)-mean(IF0(I)))/Grheobase,'b','LineWidth',1.5)
plot(time(I),(Ie0(I)+Ii0(I)+IF0(I)-mean(Ie0(I)+Ii0(I)+IF0(I)))/Grheobase,'k','LineWidth',1.5)


