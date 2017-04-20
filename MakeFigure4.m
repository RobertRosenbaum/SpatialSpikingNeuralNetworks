
clear
close all
load NetworkSimForFigure4

figure

%%%%%%%%%%%%%%%%%%%%%
% Make raster plot
%%%%%%%%%%%%%%%%%%%%%

% Size of circles plotted
mrkrsz=5;

% Time window over which to plot in ms
Tmax=2000;
Tmin=1000;

% Size of square and number of cells to plot
nplot1=20;
nplot=nplot1^2;

% indices to plot
minI=round(Ne1/2)-ceil(nplot1/2);
maxI=round(Ne1/2)+floor(nplot1/2)-1;
Iplot=find(s(1,:)>=Tmin & s(1,:)<=Tmax & s(2,:)>=minI & s(2,:)<=maxI & s(3,:)>=minI & s(3,:)<=maxI);
neuroninds=(s(2,Iplot)-minI)*nplot1+s(3,Iplot)-minI+1;
temp=randperm(nplot);
neuroninds=temp(neuroninds);

subplot(2,3,2)
plot((s(1,Iplot))/1000,neuroninds,'k.','MarkerSize',mrkrsz)


%%%%%%%%%%%%%%%%%%%%%
% Compute and plot spike count corrs
%%%%%%%%%%%%%%%%%%%%%

% counting window size
winsize=250;

% bins for distances
distbins=[0 .025:.025:.45];

% number of cells to sample
nc=5000;

% Edges for histogram
edgest=0:winsize:T;
edgesi=(1:Ne+1)-.01;
edges={edgest,edgesi};

% Find excitatory spikes,
% store into s0, which has the structure
% needed for hist3
Is=find(s(2,:)>0);
s0=zeros(numel(Is),2);
s0(:,1)=s(1,Is);
s0(:,2)=(s(2,Is)-1)*Ne1+s(3,Is);

% Get 2D histogram of spike indices and times
counts=hist3(s0,'Edges',edges);

% Get rid of edges, 
% the last element in each
% direction is zero
counts=counts(ceil(Tburn/winsize):end-1,1:end-1);

% Find neurons with rates >=1Hz
Igood=find(mean(counts)/winsize>1/1000);

% Randomly choose nc cells to use
% their indices are stored in Inds
temp=randperm(numel(Igood),nc);
Inds=Igood(temp);

% Find their locations in the square
xInds=floor((Inds-1)/Ne1+1);
yInds=mod(Inds-1,Ne1)+1;
xlocs=xInds./Ne1;
ylocs=yInds./Ne1;

% Store their spike counts
counts=counts(:,Inds);

% Compute their pairwise distances
[xlocs1,xlocs2]=meshgrid(xlocs,xlocs);
[ylocs1,ylocs2]=meshgrid(ylocs,ylocs);
distfun=@(x1,y1,x2,y2)(sqrt(min(abs(x1-x2),1-abs(x1-x2)).^2+min(abs(y1-y2),1-abs(y1-y2)).^2));
distances=distfun(xlocs1,ylocs1,xlocs2,ylocs2);
clear xlocs1 xlocs2 ylocs1 ylocs2;

% Compute their pairwise correlations
C=corrcoef(counts);

% Only keep the lower-half of the correlation matrix
% store into SCorree
[II,JJ]=meshgrid(1:size(C,1),1:size(C,1));
temp=find(II>JJ & isfinite(C));
SCorree=C(temp);

% Do same with distances
distances=distances(temp);

clear C II JJ temp;

% histogram distances
[~,I]=histc(distances,distbins);
% Get rid of meaningless last bin
distbins=distbins(1:end-1);

% Compute mean correlation and stderr of correlation
% over each distance bin
SCcorrSimMean=zeros(size(distbins));
SCcorrSimErr=zeros(size(distbins));
for j=1:numel(distbins)
    SCcorrSimMean(j)=mean(SCorree(I==j & SCorree>-2));
    SCcorrSimErr(j)=std(SCorree(I==j & SCorree>-2))./sqrt(nnz((I==j & SCorree>-2)));
end

% Plot correlations
subplot(2,3,3)
plot(distbins,SCcorrSimMean,'k','Linewidth',1.5)
hold on
axis([-.02 .47 min(SCcorrSimMean-SCcorrSimErr)-.05e-1 max(SCcorrSimMean+SCcorrSimErr)+.05e-1])

% Plot correlation distribution
[h,b]=hist(SCorree,100);
h=h./trapz(b,h);
subplot(2,3,5)
plot(b,h,'k','Linewidth',2)
temp=axis;
temp(1)=-1;
temp(2)=1;
axis(temp);


%%%%%%%%%%
% Plot current covariance and current distributions
%%%%%%%%%%

% Function to compute distances
distfun=@(x1,y1,x2,y2)(sqrt(min(abs(x1-x2),1-abs(x1-x2)).^2+min(abs(y1-y2),1-abs(y1-y2)).^2));

% distance bins
bd=[.025:.025:.708];


% Time constant of filter
tauKc=15;

% Low-pass filter
Kc=exp(-abs(-5*tauKc:dt:5*tauKc)/tauKc);
Kc=Kc/sum(Kc);

% Two-dim filter
Kcc=zeros(size(IF,1)*2+1,numel(Kc));
Kcc(size(IF,1)+1,:)=Kc;

% Low-pass filter currents
IF0=conv2(IF,Kcc,'same');
Ie0=conv2(Ie,Kcc,'same');
Ii0=conv2(Ii,Kcc,'same');
clear Kcc;

% Get rid of beginning and end, 
% which are corrupted by initial transient
% and by filtering
IF0=IF0(:,round(.5*Tburn/dt):end-round(.5*Tburn/dt));
Ie0=Ie0(:,round(.5*Tburn/dt):end-round(.5*Tburn/dt));
Ii0=Ii0(:,round(.5*Tburn/dt):end-round(.5*Tburn/dt));

% locations of recorded neurons
xlocs=Irecord(1,:)/Ne1;
ylocs=Irecord(2,:)/Ne1;

% Distances between neurons
[x1,x2]=meshgrid(xlocs,xlocs);
[y1,y2]=meshgrid(ylocs,ylocs);
distances=distfun(x1(:),y1(:),x2(:),y2(:));

% Compute covariance matrix between
% all ffwd and all rec inputs
AllCovs=cov([IF0; Ie0+Ii0]');

% Get ffwd-rec covariances
FRCovs=AllCovs(1:nrecord0,(nrecord0+1):end);
FRCovs=FRCovs(:);

% Get ffwd-ffwd covs
FFCovs=AllCovs(1:nrecord0,1:nrecord0);
FFCovs=FFCovs(:);

% Get rec-rec covs
RRCovs=AllCovs(nrecord0+1:end,nrecord0+1:end);
RRCovs=RRCovs(:);

TotalCovs=FFCovs+RRCovs+2*FRCovs;

% Compute mean and stderr of all covariances
% over each distance bin
[~,I]=histc(distances,bd);
mFRCovs=zeros(numel(bd)-1,1);
errFRCovs=zeros(numel(bd)-1,1);
mFFCovs=zeros(numel(bd)-1,1);
errFFCovs=zeros(numel(bd)-1,1);
mRRCovs=zeros(numel(bd)-1,1);
errRRCovs=zeros(numel(bd)-1,1);
mTotalCovs=zeros(numel(bd)-1,1);
errTotalCovs=zeros(numel(bd)-1,1);
for j=1:numel(bd)-1
    mFRCovs(j)=mean(FRCovs(I==j));
    errFRCovs(j)=std(FRCovs(I==j))./sqrt(nnz((I==j)));

    mFFCovs(j)=mean(FFCovs(I==j));
    errFFCovs(j)=std(FFCovs(I==j))./sqrt(nnz((I==j)));

    mRRCovs(j)=mean(RRCovs(I==j));
    errRRCovs(j)=std(RRCovs(I==j))./sqrt(nnz((I==j)));

    mTotalCovs(j)=mean(TotalCovs(I==j));
    errTotalCovs(j)=std(TotalCovs(I==j))./sqrt(nnz((I==j)));
end
bd=bd(1:end-1);

% Plot covariances
subplot(2,3,4)
plot(bd,mFRCovs/mFFCovs(1),'m','LineWidth',2)
hold on
plot(bd,mFFCovs./mFFCovs(1),'b','LineWidth',2)
plot(bd,mRRCovs./mFFCovs(1),'r','LineWidth',2)
plot(bd,mTotalCovs./mFFCovs(1),'k','LineWidth',2)
axis tight


% Current distributions
Grheobase=CalcRheoBaseEIF(Cm(1),gl(1),vl(1),DeltaT(1),vT(1),vth(1),vl(1),0,10,1000,dt,20);
subplot(2,3,6)
[h,b]=hist((Ie(:)+IF(:))/Grheobase,100);
h=h./trapz(b,h);
plot(b,h,'b','LineWidth',2)
hold on
plot(mean(Ie(:)+IF(:))/Grheobase+[0 0],[0 .02],'b','LineWidth',2)
[h,b]=hist(Ii(:)/Grheobase,100);
h=h./trapz(b,h);
plot(b,h,'r','LineWidth',2)
plot(mean(Ii(:))/Grheobase+[0 0],[0 .02],'r','LineWidth',2)
[h,b]=hist((Ie(:)+IF(:)+Ii(:))/Grheobase,100);
h=h./trapz(b,h);
plot(b,h,'k','LineWidth',2)
plot(mean(Ii(:)+Ie(:)+IF(:))/Grheobase+[0 0],[0 .02],'k','LineWidth',2)
axis tight
temp=axis;
temp(1)=-29;
temp(2)=19;
axis(temp)


