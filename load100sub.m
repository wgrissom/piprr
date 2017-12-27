function [goodFeatures,maps,fnamedirs,noiseFeatures,indices,featSNR] = load100sub(indices,nHood,sqcube,signf,nFSDs,ACCRE,addnoise,SNR,skipcoils,crossfeats,dukeorella,genderf)
%load data from 100 subs from Oct 2016 simulations.
% assumes that # of subs is however many it finds in the folder and that
% they all have a predetermined number of slices
%same as loadData_10JDI but also allows products of features
%other features.
%nFSDs = # of Fourier shape descriptors to include per slice

noiseFeatures=[];
if ~exist('dukeorella','var')
    dukeorella=[];
end

if ~exist('addnoise','var')
    addnoise=0;
end
if ~exist('skipcoils','var')
    skipcoils=0; %don't skip any coil features
elseif skipcoils==2
    keepcoils=col([1:2:8; 10:2:16; 17:2:24]);
elseif skipcoils==3
    keepcoils=col([1:3:8 10:3:16 19:3:24]);
elseif skipcoils==4
    keepcoils=col([1:4:8; 11:4:16; 17:4:24]);
elseif skipcoils==6
    keepcoils=col([1;7;13;19]);
elseif skipcoils==8
    keepcoils=[1;12;22];
elseif skipcoils==12
    keepcoils=[1;21];
elseif skipcoils==23
    keepcoils=12;
elseif skipcoils==24
    keepcoils=[];
end

% get the set of 'acceptable' shims; associate with their slice's features
%nHood  % B1+ DFT neighborhood
Nc = 24; % # tx coils
% this time.. don't run msShim to start.
maps = {};
if ACCRE
    if isempty(dukeorella) || logical(strcmp(dukeorella,'both'))
        fnamedirs = dir('./sim50Sub24ch/*24ch_10gVOP.mat');
    elseif regexp(dukeorella,'Duke')
        fnamedirs = dir('./sim50Sub24ch/Duke*24ch_10gVOP.mat');
    elseif regexp(dukeorella,'Ella')
        fnamedirs = dir('./sim50Sub24ch/Ella*24ch_10gVOP.mat');
    end
else
    if isempty(dukeorella) || logical(strcmp(dukeorella,'both'))
        fnamedirs = dir('~/sim50Sub24ch/*24ch_10gVOP.mat');
    elseif regexp(dukeorella,'Duke')
        fnamedirs = dir('~/sim50Sub24ch/Duke*24ch_10gVOP.mat');
    elseif regexp(dukeorella,'Ella')
        fnamedirs = dir('~/sim50Sub24ch/Ella*24ch_10gVOP.mat');
    end
end
% if strcmp(dukeorella,'Duke')
%     fnamedirs=fnamedirs(1:end/2);
% elseif strcmp(dukeorella,'Ella')
%     fnamedirs=fnamedirs(end/2+1:end);
% end
disp 'loading All Solutions and Features...'
goodFeatures = [];
slpersub = 31; %z slices per subject
slicetosub=1:length(fnamedirs);
slicetosub=sort(col(repmat(slicetosub,[slpersub,1])));
for subind = 1:length(fnamedirs)
    if any(slicetosub(indices)==subind) %only load sub if it's part of this
        % load this b1 case
        namebase = fnamedirs(subind).name;
        namebase = namebase(1:end-11);
        gender=isempty(regexp(namebase,'Duke'));
        if ACCRE
            x = load(['./sim50Sub24ch/' namebase '.mat']);
        else
            x = load(['~/sim50Sub24ch/' namebase '.mat']);
        end
        x.maps.mask = x.maskROI_out;
        x.dimxyz = size(x.maskROI_out);

%         %normalize the B1 maps:
%         sos_b1=sqrt(sum(abs(x.B1p_out.^2),4));
%         for sl=1:size(x.B1p_out,3)
%             sos_b1sl=sos_b1(:,:,sl);
%             x.B1p_out(:,:,sl,:)=x.B1p_out(:,:,sl,:)./mean(sos_b1sl(x.maskROI_out(:,:,sl)));
%         end
        
        x.maps.b1 = x.B1p_out;
        if addnoise
            %add in noise that gets transferred to noisy features, but save a
            %copy of clean b1 map (and clean features)
            for aa=1:size(x.B1p_out,3) %loop over slices.
                maxsig=squeeze(max(x.B1p_out(:,:,aa,:),[],4)); %peak b1 over coils for this slice
                meansig=mean(abs(maxsig(abs(maxsig)>.1*max(abs(maxsig(:)))))); %mean b1 signal for this slice               
                x.maps.noisyb1(:,:,aa,:)=x.B1p_out(:,:,aa,:)+(meansig./SNR).*(1./sqrt(2)).*(randn(size(x.B1p_out(:,:,aa,:)))+1i.*(randn(size(x.B1p_out(:,:,aa,:)))));
            end
        end
        maps{subind} = x.maps;
        if ACCRE
            load(['./sim50Sub24ch/' namebase '_10gVOP.mat'],'Sv','nSv');
        else
            load(['~/sim50Sub24ch/' namebase '_10gVOP.mat'],'Sv','nSv');
        end
        %% get SAR
        % get the average SAR matrix, and append it to the rest
        Sv(:,:,end+1) = Sv(:,:,1)*nSv(1);
        for ii = 2:size(Sv,3)-1
            Sv(:,:,end) = Sv(:,:,end) + Sv(:,:,ii)*nSv(ii);
        end
        Sv(:,:,end) = Sv(:,:,end)/sum(nSv);
        % scale it up 5x since the max local can be e.g. 20 W/kg but the max ave can only be 4 W/kg
        Sv(:,:,end) = 5*Sv(:,:,end);
        C = buildSARReg(sum(Sv,3),1);
        maps{subind}.C = C./(size(C,1)/Nc); %tack on the SAR vops for this subject, normalize by numVOPS
        maps{subind}.Sv =Sv;
        x.Nsl=size(x.maskROI_out,3);
        dims = size(x.maskROI_out);
        x.dimxy=dims(1:2);
        for ii = 1:x.Nsl
            if any(indices == (subind-1)*x.Nsl + ii) % include this slice if its on the list
                % first, check if the slice has enough voxels to be
                % considered...
                if sum(col(x.maskROI_out(:,:,ii)))>=20 %discard if less than 20 voxels

                % extract this slice's features
                % include central DFT coeffs for each coil, normalized to coil 1's DC
                % coefficient
                sliceFeatures = ii; % store slice position
                %sliceFeatures = x.notEmptyInds(ii)*([streq(orientations{ori},'axial');...
                %    streq(orientations{ori},'sagittal');...
                %    streq(orientations{ori},'coronal')]);
                % calculate the 'width' of the mask in x and y & z
                [xi,yi,zi] = ndgrid(1:x.dimxyz(1),1:x.dimxyz(2),1:x.dimxyz(3));
                
                %***TESTING*** JDI
                if exist('genderf') && genderf
                sliceFeatures = [sliceFeatures; gender];
                end
                %global brain dimensions:
%                 sliceFeatures = [sliceFeatures;std(xi(x.maskROI_out));std(yi(x.maskROI_out));std(zi(x.maskROI_out))];
%                 % calculate centroid of mask
%                 sliceFeatures = [sliceFeatures;mean(xi(x.maskROI_out));mean(yi(x.maskROI_out));mean(zi(x.maskROI_out))];
                
                % slice-spef dimensions:
                sliceFeatures = [sliceFeatures;std(xi(x.maskROI_out(:,:,ii)));std(yi(x.maskROI_out(:,:,ii)))];
                sliceFeatures = [sliceFeatures;mean(xi(x.maskROI_out(:,:,ii)));mean(yi(x.maskROI_out(:,:,ii)))];
                
%                 [xi,yi] = meshgrid(1:x.dimxy(1),1:x.dimxy(2));
%                 sliceFeatures = [sliceFeatures;std(xi(x.maskROI_out(:,:,ii)));std(yi(x.maskROI_out(:,:,ii)))];
%                 
%                 % calculate centroid of mask
%                 sliceFeatures = [sliceFeatures;mean(xi(x.maskROI_out(:,:,ii)));mean(yi(x.maskROI_out(:,:,ii)))];
                %% add Fourier shape descriptors
                if nFSDs>0
                    [~,fdes]=fourierShapeDescriptor(x.maskROI_out(:,:,ii),nFSDs);
                    sliceFeatures=[sliceFeatures; fdes(:)];
                end
                %% get neighborhood B1 coeffs
                
                if nHood > 0
                    for jj = 1:Nc
                        %only add this coil's features if we're not dropping
                        %coils OR this is a non-dropped coil
                        if ~skipcoils
                            cont=1;
                        elseif any(keepcoils==jj)
                            cont=1;
                        else
                            cont=0;
                        end
                        if cont
                            if jj==1
                                nonb1feats=sliceFeatures;
                                noise=[];
                                tmp3=[];
                            end
                            %non-noisy b1 features
                            tmp = x.maps.b1(:,:,ii,jj);
                            tmp = fftshift(fft2(fftshift(tmp)));
                            tmp = col(tmp(ceil(x.dimxy(1)/2)+1-(nHood-1)/2:ceil(x.dimxy(1)/2)+1+(nHood-1)/2,...
                                ceil(x.dimxy(2)/2)+1-(nHood-1)/2:ceil(x.dimxy(2)/2)+1+(nHood-1)/2));
                            sliceFeatures = [sliceFeatures; tmp];

                            %noisy b1 features
                            if addnoise
                                tmp2 = x.maps.noisyb1(:,:,ii,jj);
                                tmp2 = fftshift(fft2(fftshift(tmp2)));
                                tmp2 = col(tmp2(ceil(x.dimxy(1)/2)+1-(nHood-1)/2:ceil(x.dimxy(1)/2)+1+(nHood-1)/2,...
                                    ceil(x.dimxy(2)/2)+1-(nHood-1)/2:ceil(x.dimxy(2)/2)+1+(nHood-1)/2));
                                tmp3 = [tmp3; tmp2];
                                %calc feature SNR resulting from noise:
                                noise=[noise;tmp2-tmp]; 
                                if jj==Nc % if last coil add the noisy data
                                    featSNR=mean(abs(tmp3./noise));
                                    %add noise of featSNR to remaining (non-b1)
                                    %features.
                                    noisyfeats = [nonb1feats + (nonb1feats./featSNR).*(1./sqrt(2)).*(randn(size(nonb1feats))+1i.*randn(size(nonb1feats))); tmp3];
                                end
                            end
                        end
                    end
                end
                if ~exist('sqcube') || ~sqcube
                    if signf
                        sliceFeatures = [sliceFeatures; sign(sliceFeatures)];
                    end
                else %use squared and cubed slice features
                    if nHood>0
                        sliceFeatures = [sliceFeatures; sign(sliceFeatures(end-2*Nc:end)); sliceFeatures.^2; sliceFeatures.^3];
                    else
                        sliceFeatures = [sliceFeatures; sign(sliceFeatures); sliceFeatures.^2; sliceFeatures.^3];
                    end
                end
                % ditch the sign features or we'll have too many... but
                % introduce cross feats:
                if crossfeats
                    crossFeats=sliceFeatures*sliceFeatures.';
                    if addnoise
                        noisycrossfeats=noisyfeats*noisyfeats.';
                    end
                    %ditch the lower half of triangular matrix
                    crossMask=logical(triu(ones(size(crossFeats))));
                    crossFeats=crossFeats(crossMask);
                    if addnoise
                        noisycrossfeats=noisycrossfeats(crossMask);
                    end
                    sliceFeatures=[sliceFeatures; crossFeats];
                    if addnoise
                        noisyfeats=[noisyfeats; noisycrossfeats];
                    end
                end
%                 if addnoise
%                     nonoisefeats = [nonoisefeats sliceFeatures];
%                     sliceFeatures = sliceFeatures + (sliceFeatures./SNR).*randn(size(sliceFeatures));
%                 end
                goodFeatures = [goodFeatures sliceFeatures];
                if addnoise
                    noiseFeatures =[noiseFeatures noisyfeats];
                end
                % append the final set of features to the running
                % list
                else %remove this index from the list
                    indices(find(indices==(subind-1)*x.Nsl + ii))=[];
                end
            end
        end
    end
end
nFeatures = size(goodFeatures,1);
Nsl=size(goodFeatures,2);
fprintf('Loaded a total of %d features and %d slices\n',nFeatures,Nsl);
fnamedirs=fnamedirs.name;