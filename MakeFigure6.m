
clear
close all

%%%%% First for L23
load NetworkSimForFigure6L4toL23.mat

rng(1)

figure


%%%%%%%%%%%%%%%%%%%%%
% Compute and plot spike count corrs
%%%%%%%%%%%%%%%%%%%%%

% counting window size
winsize=250;

% bins for distances
distbins=[0 .05:.05:.6 2];

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
subplot(1,2,2)
plot(distbins,SCcorrSimMean,'k','Linewidth',1.5)
hold on

% Plot correlation distribution
[h,b]=hist(SCorree,100);
h=h./trapz(b,h);
subplot(1,2,1)
plot(b,h,'k','Linewidth',2)
hold on
temp=axis;
temp(1)=-1;
temp(2)=1;
axis(temp);


%%%%% Now for L4
clear s;
load NetworkSimForFigure6LGNtoL4.mat

s=s0;
clear s0;

%%%%%%%%%%%%%%%%%%%%%
% Compute and plot spike count corrs
%%%%%%%%%%%%%%%%%%%%%

% counting window size
winsize=250;

% bins for distances
distbins=[0 .05:.05:.6 2];

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
subplot(1,2,2)
plot(distbins,SCcorrSimMean,'Color',[.5 .5 .5],'Linewidth',1.5)

% Plot correlation distribution
[h,b]=hist(SCorree,100);
h=h./trapz(b,h);
subplot(1,2,1)
plot(b,h,'Color',[.5 .5 .5],'Linewidth',2)
axis tight;
temp=axis;
temp(1)=-1;
temp(2)=1;
axis(temp);

