
clear
close all
load NetworkSimForFigure4;

rng(1)

%%%%%%%%
%%%%%%%%

% Proportion of neurons in population F
qF=NF/N;

% Window size over which to compute correlations
winsize=250;

% Fourier transforum of PSC kernels
taujitter=0;
etaF=@(f)((1./(1+2*pi*sqrt(-1)*f*tausynF)).*exp(-4*pi^2*f.^2*taujitter^2));
etae=@(f)(1./(1+2*pi*sqrt(-1)*f*tausyne));
etai=@(f)(1./(1+2*pi*sqrt(-1)*f*tausyni));

% Mean-field connectivity at zero spatial Fourier mode
% as a funtion of temporal frequency (see Supp Materials)
weebar=@(f)(wee0*etae(f));
weibar=@(f)(-wei0*etai(f));
wiebar=@(f)(wie0*etae(f));
wiibar=@(f)(-wii0*etai(f));
weFbar=@(f)(sqrt(qF)*jeF*peF0*etaF(f));
wiFbar=@(f)(sqrt(qF)*jiF*piF0*etaF(f));

% Matrix-form of connectivity and ffwd input in Fourier domains
FFbar=@(f)((rF*[abs(weFbar(f)).^2 weFbar(f).*conj(wiFbar(f)); wiFbar(f).*conj(weFbar(f)) abs(wiFbar(f)).^2]));
Wbar=@(f)([weebar(f) weibar(f); wiebar(f) wiibar(f)]);
FF=@(n1,n2,f)((rF*[abs(weFbar(f)).^2 weFbar(f).*conj(wiFbar(f)); wiFbar(f).*conj(weFbar(f)) abs(wiFbar(f)).^2])*exp(-2*pi^2*(n1.^2+n2.^2)*2*sigmaeF^2));
W=@(n1,n2,f)([weebar(f)*exp(-2*pi^2*(n1^2+n2^2)*sigmaee^2) weibar(f)*exp(-2*pi^2*(n1^2+n2^2)*sigmaei^2); wiebar(f)*exp(-2*pi^2*(n1^2+n2^2)*sigmaie^2) wiibar(f)*exp(-2*pi^2*(n1^2+n2^2)*sigmaii^2)]);

% Inverse transfer fcn, 1/gain on the diagonal
% Gains are computed by fitting rates and inputs
% to a threshold quadratic, then funding derivative
% at the mean rate (see Supp Materials).
ge=.0122;
gi=.0174;
Linv=@(n1,n2)(([1/ge 0; 0 1/gi]));

% Spike-spike cross-spectral matrix in spatial Fourier domain
% Directly fro Eqn S.49 in Supp Materials.
SS=@(n1,n2,f)((((inv(Linv(n1,n2)-sqrt(N)*W(n1,n2,f))*FF(n1,n2,f)*inv(Linv(n1,n2)-sqrt(N)*W(n1,n2,f)')))));


% Uniformly sample positions in space, to compute all
% pairwise distances and pairwise correlations.
U1=0:.002:1;
U2=0:.002:1;
dists=sqrt(min(U1(:),1-U1(:)).^2+min(U2(:),1-U2(:)).^2);

% For all pairs of spatial freqs up to maxfreq
maxfreq=20;
Cee=zeros(size(U1));
for n1=-maxfreq:maxfreq
for n2=-maxfreq:maxfreq    
    
    % Get cross-spectral matrix at zero temporal freq
    temp=SS(n1,n2,0);
    
    % The spike count covariance over large windows 
    % between excitatory neurons is computed by summing the 
    % two-dimensional Fourier series of the EE (i.e., (1,1)) 
    % component of the cross-spectral matrix.
    % This line adds the n1,n2 term of this sum.
    Cee=Cee+(temp(1,1)*exp(2*pi*sqrt(-1)*(n1*U1+n2*U2)));
    
end
end

% Get distances up to .7
% note that max distance is
% sqrt(1/2)=.7071 because 
% boundaries are periodic
I=find(dists>.7,1);
dists=dists(1:I);
Cee=Cee(1:I);
DistancesTh=dists;

% Count covariance over large windows is approximately Cee*winsize
% Mean count is rate*winsize.
% Count correlation is approx. count covariance over mean count
% if spike trains are approx. Poisson
CorrsTh=real(Cee)*winsize/(winsize*reSim/1000);
save TheoreticalCorrsForFig4.mat DistancesTh CorrsTh;

