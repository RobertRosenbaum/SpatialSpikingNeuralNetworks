
clear
close all

path(path,'./gpfa-master')

load NetworkSimForFigure6L4toL23.mat

rng(1)

figure

% Start gpfa algorithm
startup()


nc=500;
dtP=50;
winsize=250;
nbins=round(winsize/dtP);
ntrials=round((T-Tburn)/winsize);


% Edges for histogram
edgest=0:dtP:T;
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
counts=counts(ceil(Tburn/dtP):end-2,1:end-1);

% Find neurons with rates >=1Hz
Igood=find(mean(counts)/dtP>1/1000);

% Randomly choose nc cells to use
% their indices are stored in Inds
temp=randperm(numel(Igood),nc);
Inds=Igood(temp);

Y=reshape(counts(:,Inds)',nc,nbins,ntrials);
model = GPFA();
model = model.fit(Y, 1, 'hist');

[R0,X]=model.residCov(Y,1);  % By trial

ResidCorrMatrix=corrcov(R0);

LatentCovMatrix=(model.C)*(model.C)';

% Find their locations in the square
xInds=floor((Inds-1)/Ne1+1);
yInds=mod(Inds-1,Ne1)+1;
xlocs=xInds./Ne1;
ylocs=yInds./Ne1;


% Compute their pairwise distances
[xlocs1,xlocs2]=meshgrid(xlocs,xlocs);
[ylocs1,ylocs2]=meshgrid(ylocs,ylocs);
distfun=@(x1,y1,x2,y2)(sqrt(min(abs(x1-x2),1-abs(x1-x2)).^2+min(abs(y1-y2),1-abs(y1-y2)).^2));
distances=distfun(xlocs1,ylocs1,xlocs2,ylocs2);
clear xlocs1 xlocs2 ylocs1 ylocs2;

[II,JJ]=meshgrid(1:nc,1:nc);

ResidCorrs=ResidCorrMatrix(II>JJ);
LatentCovs=LatentCovMatrix(II>JJ);
distances=distances(II>JJ);

% Distance on the unit square
distbins=[0 .1:.1:.5 10];%[0 .1:.1:.5 10];

[~,I]=histc(distances,distbins);
% Get rid of meaningless last bin
distbins=distbins(1:end-1);

% Compute mean correlation and stderr of correlation
% over each distance bin
ResidCorrMean=zeros(size(distbins));
ResidCorrErr=zeros(size(distbins));
LatentCovMean=zeros(size(distbins));
LatentCovErr=zeros(size(distbins));
for j=1:numel(distbins)
    ResidCorrMean(j)=mean(ResidCorrs(I==j & ResidCorrs>-2));
    ResidCorrErr(j)=std(ResidCorrs(I==j & ResidCorrs>-2))./sqrt(nnz((I==j & ResidCorrs>-2)));

    LatentCovMean(j)=mean(LatentCovs(I==j & LatentCovs>-2));
    LatentCovErr(j)=std(LatentCovs(I==j & LatentCovs>-2))./sqrt(nnz((I==j & LatentCovs>-2)));
end

errorbar(distbins*10,ResidCorrMean,ResidCorrErr,'Color',[.5 .5 .5],'LineWidth',2)

