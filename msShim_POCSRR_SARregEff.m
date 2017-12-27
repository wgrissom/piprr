% Script to do multi-slice RF shimming
% Copyright Zhipeng Cao and Will Grissom, Vanderbilt University, 2015.

% TO IMPLEMENT:
% 1) Load slice/subject indices for train cases as in compare_kNN
% 2) Get features for train cases
% 3) Normalize features for train cases; add some feature transformations
% 4) Build Feature projector for train cases, pass to POCS ms-shim
    
[dimxy(1),dimxy(2),Nsl,Nc] = size(mapsTrain.b1); % number of physical coils (Ncoils)

% normalize b1 by median so that regularization params can be
% meaningfully reported
mapsTrain.b1 = mapsTrain.b1./median(abs(mapsTrain.b1(repmat(mapsTrain.mask,[1 1 1 Nc]))));

algp.beta = 0; % RF power regularization
algp.tol = 1e-6; %1-0.9999; % stopping parameter
if ~isfield(algp,'noCG')
    algp.noCG = false;
end
algp. step =0.1;
algp.nRandStart = -1; % only do CP initialization
algp.ncgiters = ncgiters; %3; % number of cg iters for each RF update
algp.dofigs = 0; % whether to show pattern figure every 100 iters
%algp.holeThresh = 0.25; % hole detection threshold
%algp.holeReg = holeReg; % whether to monitor holes
% Set up SAR regularization
% algp.betaSAR=betaSAR; %10^-2; % SAR regularization weight
algp.cutoff=cutoff;
if exist('Fproj','var') 
    algp.Fproj = Fproj; % projector onto space of predictable solutions
end

[rfPOCS,errAllPOCS,SAR_AllPOCS,mAllPOCS,randAmpPhsPOCS,algp] = msShim_randStart_POCSRR_SARregEff(mapsTrain,algp);