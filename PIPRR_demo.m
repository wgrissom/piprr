function PIPRR_demo(kk,nHood,betaSAR,ncgiters,cutoff,maxiters,debugging,nFSD,ACCRE,regperc,crossfeats,decPR,skipcoils,dukeorella,addnoise,SNR,svdreg,genderf,iterative)
% demo script to run RF shim PIPRR

%set up parameters:
nHood = 1;      %neighborhood size of B1 map Fourier coefficients 
                    %to use as features.
betaSAR = 0.01; %controls strength of SAR regularization
ncgiters = 3;   %number of CG iterations used in MLS shim design
nFSD = 3;       %number of Fourier shape descriptors used as features
regperc = 1;    %controls the strength of KRR regularization

%Nccomp = number of coils to compress to
%nHood =  fourier neighborhood over which to extract B1+ coeff features%POCSreg = value of POCs regularization
%sqcube = whether to include cube and square of features
%subout = 1 if we are leaving an entire subject out for the testing.
%decPR = 1 if we're decreasing POCSreg as we go...
%set kk= whatever (K-)fold we asked for.
% this is WITH 1st order crossfeats.
% also validation-based cutoff of POCS iters.
% algp.bSARperc=bSARperc;
% iterative = whether to do the iterative ML procedure (that works) or the
% non-iterative straight linear ridge reg. (fatally flawed)

if ~exist('iterative','var')
    iterative=1;
elseif iterative==0
    algp.ncgiters=20; %up this so that we can get good/realistic sols to training data
    maxiters=1;
end
if ~exist('dukeorella','var') 
    dukeorella=[];
end
    if ~exist('genderf','var')
    genderf=0;
end
algp.betaSAR=betaSAR;
algp.bSARperc=betaSAR;
starttime=tic
if exist('decPR','var')
    algp.decPR=decPR;
else
    algp.decPR=0;
end
%*NOTE* k-folds are different in this script bc it operates over subjects,
%not single slices!
if ~exist('debugging','var')
    debugging=0;
end
algp.maxIters = maxiters; % max iters before exiting
if ACCRE
cd /home/iannijd/EPI_TrACR/fessler/irt/
setup
cd /home/iannijd/RFRF/SARsim/
else
addpath(genpath('~/Copy/research/RFdictionaries/tensub/'));
addpath(genpath('~/Copy/research/tensorflow'));
addpath(genpath('~/Copy/research/misc_functions/'));
rmpath(genpath('~/Copy/research/misc_functions/SheffieldML-multigp-c5c1d8e/'));
addpath(genpath('~/Copy/research/RFdictionaries/SARreg/'))
addpath(genpath('~/Copy/acptx/ktpoints'))
end
signf=0;
sqcube=0;
K = 10; % K-fold cross validation
Nsl = 31; % # slices in each subject's maps
if ~isempty(dukeorella)
    if regexp(dukeorella,'both')
        if ACCRE
            Nsubj=length(dir('./sim50Sub24ch/*24ch_10gVOP.mat'));
        else
            Nsubj=length(dir('~/Copy/sim50Sub24ch/*24ch_10gVOP.mat'));
        end
    else
        Nsubj=50;
    end
else
    Nsubj = 100; % # subjects
end
%filename and name of 1st k-fold file...
algp.fstring=['trePRR100sub_SARreg' num2str(betaSAR) '_kfold' num2str(kk) 'Hood' num2str(nHood) 'ncg' num2str(ncgiters) 'cut' num2str(cutoff) 'nFSD' num2str(nFSD) 'POCSreg' num2str(regperc) 'cross' num2str(crossfeats) 'decPR' num2str(decPR) 'dprop' num2str(exist('decPRprop')) 'skip' num2str(skipcoils) dukeorella num2str(Nsubj) 'SNR' num2str(SNR) 'gen' num2str(genderf) 'it' num2str(iterative) '.mat'];
algp.fstringItr=['trePRR100sub_SARreg' num2str(betaSAR) '_kfold1Hood' num2str(nHood) 'ncg' num2str(ncgiters) 'cut' num2str(cutoff) 'nFSD' num2str(nFSD) 'POCSreg' num2str(regperc) 'cross' num2str(crossfeats) 'decPR' num2str(decPR) 'dprop' num2str(exist('decPRprop')) 'skip' num2str(skipcoils) dukeorella num2str(Nsubj) 'SNR' num2str(SNR) 'gen' num2str(genderf) 'it' num2str(iterative) '.mat'];

Nc = 24; %* note not 36 !

% POCSreg = regularization on POCS pseudoinverse

if Nsubj==100;
    if exist('kfoldhowmuchtrain.mat')
        load('kfoldhowmuchtrain.mat');
    else
        subjkfoldInds = crossvalind('Kfold',Nsubj,Nsubj); % what k-fold each subj belongs to
        slicetosub=1:Nsubj;
        slicetosub=sort(col(repmat(slicetosub,[Nsl,1])));
        kfoldInds = subjkfoldInds(slicetosub); % which k-fold each slice belongs to.
        save('kfoldhowmuchtrain.mat','subjkfoldInds','slicetosub','kfoldInds');
    end
elseif Nsubj==50;
    if exist('kfoldhowmuchtrain50.mat')
        load('kfoldhowmuchtrain50.mat');
    else
        subjkfoldInds = crossvalind('Kfold',Nsubj,Nsubj); % what k-fold each subj belongs to
        slicetosub=1:Nsubj;
        slicetosub=sort(col(repmat(slicetosub,[Nsl,1])));
        kfoldInds = subjkfoldInds(slicetosub); % which k-fold each slice belongs to.
        save('kfoldhowmuchtrain50.mat','subjkfoldInds','slicetosub','kfoldInds');
    end
else
        if exist(['kfoldhowmuchtrain' num2str(Nsubj) '.mat'])
            load(['kfoldhowmuchtrain' num2str(Nsubj) '.mat']);
        else
        subjkfoldInds = crossvalind('Kfold',Nsubj,Nsubj); % what k-fold each subj belongs to
        slicetosub=1:Nsubj;
        slicetosub=sort(col(repmat(slicetosub,[Nsl,1])));
        kfoldInds = subjkfoldInds(slicetosub); % which k-fold each slice belongs to.
        save(['kfoldhowmuchtrain' num2str(Nsubj) '.mat'],'subjkfoldInds','slicetosub','kfoldInds');
        end
end


% init solution + error arrays
mPOCS = zeros(46,55,Nsubj*Nsl);
percErrmPOCS = zeros(length(kfoldInds),1);
solErrPOCS = zeros(length(kfoldInds),1);
mTrainPOCSRR=[];
percErrmTrainPOCSRR=[];
kktest=(kk-1)*floor(Nsubj/10)+1:kk*(floor(Nsubj/10));
kktrain=1:Nsubj;
kktrain=kktrain(~ismember(kktrain,kktest));
if ~debugging
    testInds=[];
    for ii=1:length(kktest)
        testInds = [testInds; find((kktest(ii)==kfoldInds))];
    end
    trainInds=[];
    for ii=1:length(kktrain)
        trainInds = [trainInds; find((kktrain(ii)==kfoldInds))];
        subjtrainInds = unique(slicetosub(trainInds));
    end
    nTest = length(testInds);
    nTrain = length(trainInds);
else
    testInds=[1:12];
    trainInds=[13:24];
    subjtrainInds=ones(numel(trainInds),1);
    nTest=numel(testInds);
    nTrain=numel(trainInds);
end
testInds=sort(testInds);
trainInds=sort(trainInds);
% addnoise=0;
% SNR=Inf;
if ~exist('skipcoils','var')
skipcoils=0;
end

    % load the shims and features of the slices we don't train on -
    % REDFLAG: Assumes testInds is sorted!
    [testFeatures,maps,fnamedirs,nonoisefeatstest,testInds] = load100sub(testInds,nHood,sqcube,signf,nFSD,ACCRE,addnoise,SNR,skipcoils,crossfeats,dukeorella,genderf);
    
    % load the training shims and features - REDFLAG: Assumes trainInds is
    % sorted!
    
    [trainFeatures,trainmaps,~,nonoisefeatstrain,trainInds] = load100sub(trainInds,nHood,sqcube,signf,nFSD,ACCRE,addnoise,SNR,skipcoils,crossfeats,dukeorella,genderf);

%reset train & test #s since we've discarded slices
nTrain=length(trainInds);
nTest=length(testInds);
    
% fill in subject maps missing (for all subs not in test set)
for tr_indind = 1:length(subjtrainInds)
    tr_ind = subjtrainInds(tr_indind);
    maps{tr_ind}=trainmaps{tr_ind};
end
algp.testInds=testInds;
algp.trainInds=trainInds;
%% run the POCS-RR-shim design
% first get all the maps for the design into one array
mapsTrain.b1 = [];
mapsTrain.mask = [];
for ii = 1:nTrain
    subjInd = floor(trainInds(ii)/Nsl)+1;
    slInd = rem(trainInds(ii),Nsl);
    if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end
    mapsTrain.b1 = cat(3,mapsTrain.b1,maps{subjInd}.b1(:,:,slInd,:));
    mapsTrain.mask = logical(cat(3,mapsTrain.mask,maps{subjInd}.mask(:,:,slInd)));
    mapsTrain.C{ii} = maps{subjInd}.C; %tack on the SAR vops for this slice's subject,
        % already normalized by numVOPS
end
% then get all the features for the design
% get normalized training features for nearest-neighbors methods


featMean = mean(trainFeatures,2);
featStd = std(trainFeatures,0,2);
trainFeatures = (trainFeatures-repmat(featMean,[1 nTrain]))...
    ./repmat(featStd,[1 nTrain]);
trainFeatures(isinf(trainFeatures) | isnan(trainFeatures)) = 0; % get divide by zero for some coordinates
% get normalized test features
%featMean = mean(testFeatures,2);
%featStd = std(testFeatures,0,2);
testFeatures = (testFeatures-repmat(featMean,[1 nTest]))...
    ./repmat(featStd,[1 nTest]);
testFeatures(isinf(testFeatures) | isnan(testFeatures)) = 0; % get divide by zero for some coordinates

% add a 1 to each slice's feature vector to account for offset
trainFeatures = [trainFeatures;ones(1,nTrain)];
testFeatures = [testFeatures;ones(1,nTest)];

%(just to save):
algp.featMean=featMean; 
algp.featStd=featStd;
algp.trainFeatures=trainFeatures;
algp.testFeatures=testFeatures;

% calculate projector matrix
F = trainFeatures.';
algp.F=F;
if ~exist('svdreg') || ~svdreg
    % [~,S,~]=svd(F);
    POCSreg=regperc;%(regperc.*max(abs(diag(S))));
else
    [~,S,~]=svd(F);
    POCSreg=(regperc.*max(abs(diag(S))));
end
algp.POCSreg=POCSreg;
Fproj = F*((F'*F + POCSreg*eye(size(F,2)))\F');
Finv = (F'*F + POCSreg*eye(size(F,2)))\F';

algp.compress=0;
algp.Finv=Finv;
algp.kk=kk;
if kk==1
    algp.maps=maps;
    algp.testInds=testInds;
    algp.testFeaturesPOCSshim=testFeatures;
    algp.Finv=Finv;
end

Nslpersub=Nsl;
algp.compress=0;
msShim_POCSRR_SARregEff;
% get prediction matrix back and apply to get the test RF
rfPOCS = sqz(rfPOCS);
% get A matrix
APOCS = [];
for ii = 1:Nc
%     APOCS = [APOCS; (Finv*real(rfPOCS(ii,:)).').'; (Finv*imag(rfPOCS(ii,:)).').'];
    APOCS = [APOCS; (Finv*rfPOCS(ii,:).').'];
end

%% get the training sols:
% get the errors of the POCS-RR-designed training shims:
for ii = 1:nTrain
    
    subjInd = floor(trainInds(ii)/Nslpersub)+1;
    slInd = rem(trainInds(ii),Nslpersub);
    if slInd == 0;slInd = Nslpersub;subjInd = subjInd - 1;end;
    
    tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
    for jj = 1:Nc
        tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*rfPOCS(jj,ii);
    end
    mTrainPOCSRR(:,:,end+1) = tmp;
    percErrmTrainPOCSRR(end+1) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
end
%% get the test solutions
for ii = 1:nTest
    Nslpersub=31;
    predPOCS = APOCS*testFeatures(:,ii);
    % get subject and slice indices
    subjInd = floor(testInds(ii)/Nslpersub)+1;
    slInd = rem(testInds(ii),Nslpersub);
    if slInd == 0;slInd = Nslpersub;subjInd = subjInd - 1;end;
    
    tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
    for jj = 1:Nc
        tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predPOCS(jj);
    end
    mPOCS(:,:,testInds(ii)) = tmp;
    percErrmPOCS(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
    
    % get error between best solutions and predicted solutions
%     solErrPOCS(testInds(ii)) = norm(testSols(:,ii) - predPOCS)/norm(testSols(:,ii));
end
if isfield(algp,'maps')
    algp=rmfield(algp,'maps');
end
save(algp.fstring,'APOCS','rfPOCS','percErrmTrainPOCSRR','percErrmPOCS','SAR_AllPOCS','algp','fnamedirs','nonoisefeatstest','nonoisefeatstrain','-v7.3');
toc(starttime)
