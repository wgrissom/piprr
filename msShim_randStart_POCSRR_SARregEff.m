function [rf,errAll,vopAll,mAll,randAmpPhs,algp] = msShim_randStart_POCSRR_SARregEff(maps,algp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta = algp.beta; % 0    % RF regularization parameter
if ~isfield('algp','cutoff')
    algp.cutoff=1;
end
algp.regSave=[];
algp.initreg=[];
% betaSAR=algp.betaSAR;
ncgiters=algp.ncgiters;
% compress = algp.compress; %1 if using array-compression, 0 if not...
%Nvops=algp.Nvops;
%C = %SAR regularization term from Sv & buildSARReg. (but normalized! by
%       # of VOPS
%for qwpls_pcg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build system matrix for each slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dimxy(1),dimxy(2),Nsl,Nc] = size(maps.b1);
% for idx=1:Nsl
%     c2=full(maps.C{idx});
%     c2=c2./std(c2(:));
%     maps.C{idx}=c2;
% end
%normalize VOPS for regularization purposes
% Cbig=Cbig./std(Cbig(:));
% for idx=1:Nsl
%     maps.C{idx}=Cbig(:,:,idx);
% end
for idx = 1:Nsl
    
    % reorg b1 maps into tall skinny matrix (row: space, col: coil)
    tmp = permute(squeeze(maps.b1(:,:,idx,:)),[3 1 2]);
    tmp = tmp(:,:).';
    % mask them
    tmp = tmp(logical(maps.mask(:,:,idx)),:);
    % remove coil 1 phase
    A{idx} = bsxfun(@times,exp(-1i*angle(tmp(:,1))),tmp);
%     [~,S,~]=svd(A{idx});
%     s=diag(S);
%     s_1=max(abs(s));
%     [~,S,~]=svd(full(maps.C{idx}));
%     s_2=max(abs(diag(S)));
    algp.betaSAR(idx)=algp.bSARperc;%.*s_1./s_2;
    if algp.noCG
        ARinv{idx} = (A{idx}'*A{idx} + algp.betaSAR*(maps.C{idx}'*maps.C{idx}))\A{idx}';
    end
    
    % init target phase pattern
    if ~isfield(maps,'phsinit')
        dphs{idx} = zeros(sum(col(maps.mask(:,:,idx))),1);
    else
        tmp = maps.phsinit(:,:,idx);
        dphs{idx} = tmp(col(maps.mask(:,:,idx)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do compressed design for increasing # of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a set of random amplitudes and phases to use in initializing each
% slice's target phase pattern
if ~isfield(algp,'randAmpPhs')
    randAmpPhs = rand(Nc,2,algp.nRandStart);
else
    randAmpPhs = algp.randAmpPhs;
end
nRandStart = algp.nRandStart + 2; % also add zero phase and direct sum phase
rf = zeros(Nc,nRandStart,Nsl);
errAll = zeros(nRandStart,Nsl);
vopAll = zeros(nRandStart,Nsl);
mAll = zeros(dimxy(1),dimxy(2),nRandStart,Nsl);
tol = algp.tol;
cost=Inf;
valcost=Inf;
for raIdx = 1:nRandStart
    
    % get initial phase patterns for each slice
    for slIdx = 1:Nsl
        
        % precalculate regularized pseudoinverse
        %Ainv{slIdx} = (A{slIdx}'*A{slIdx} + R)\A{slIdx}';
        
        % get random initial starting phase
        switch raIdx
            case 2
                dphs{slIdx} = zeros(size(A{slIdx},1),1);
            case 1
                dphs{slIdx} = angle(sum(A{slIdx},2));
            otherwise
                dphs{slIdx} = angle(A{slIdx}*(randAmpPhs(:,1,raIdx-2).*exp(1i*2*pi*randAmpPhs(:,2,raIdx-2))));
        end
        
    end
    
    % outer loop: continue until shims stop changing
    flag = 1; % to get loop started
    itr = 1; % iteration counter
    algp.holesExist = false;
    while flag
        tic
        rfOld = sqz(rf(:,raIdx,:));
        
        % update each slice's shims to minimize shim error
        if algp.noCG
            for slIdx = 1:Nsl
                rf(:,raIdx,slIdx) = (1-algp.step)*rf(:,raIdx,slIdx) + ...
                    algp.step*(ARinv{slIdx}*exp(1i*dphs{slIdx}));
            end
        else
            for slIdx = 1:Nsl
                % take a couple CG iterations
                xS = qpwls_pcg(rf(:,raIdx,slIdx),A{slIdx},1,...
                    exp(1i*dphs{slIdx}),0,sqrt(algp.betaSAR(slIdx))*maps.C{slIdx},1,ncgiters,ones(Nc,1));
                %xS = qpwls_pcg_ctc(rf(:,raIdx,slIdx),A{slIdx},1,...
                %    exp(1i*dphs{slIdx}),0,beta*(maps.R{slIdx}'*maps.R{slIdx}),1,algp.ncgiters,ones(Nc,1));
                rf(:,raIdx,slIdx) = xS(:,end);
            end
        end
        % project shims onto space of predictable shims
        if isfield(algp,'Fproj')
            if length(valcost)>=2 && rem(itr-1,10)==0 && algp.decPR %update POCSreg
                algp.regSave(end+1)=algp.POCSreg;
                algp.POCSreg=algp.POCSreg.*.99;%*((valcost(end)./cost(end)));
                algp.Fproj = algp.F*((algp.F'*algp.F + algp.POCSreg*eye(size(algp.F,2)))\algp.F');
                if algp.kk==1
                    algp.Finv = (algp.F'*algp.F + algp.POCSreg*eye(size(algp.F,2)))\algp.F';
                end
            elseif length(valcost)>=2 && rem(itr-1,10)==0 && isfield(algp,'decPRprop') %update POCSreg
                algp.regSave(end+1)=algp.POCSreg;
                algp.POCSreg=algp.POCSreg.*.99*((valcost(end)./cost(end)));
                algp.Fproj = algp.F*((algp.F'*algp.F + algp.POCSreg*eye(size(algp.F,2)))\algp.F');
                if algp.kk==1
                    algp.Finv = (algp.F'*algp.F + algp.POCSreg*eye(size(algp.F,2)))\algp.F';
                end
            end
            
            saverflast =  rf;
            
            for ii = 1:Nc
                % project real and imag parts of each slice's shim for this
                % coil onto the space of predictable shims, assuming that
                % the prediction weights are given by fitting the
                % features to this training data across slices
                rf(ii,raIdx,:) = algp.Fproj*sqz(rf(ii,raIdx,:));
            end
        end
        
        %calculate the KRR cost 
        rfPOCS = sqz(rf);
        % get A matrix
        APOCS = [];
        for ii = 1:Nc
            %         APOCS = [APOCS; (algp.Finv*real(rfPOCS(ii,:)).').'; (algp.Finv*imag(rfPOCS(ii,:)).').'];
            APOCS(ii,:) = (algp.Finv*rfPOCS(ii,:).');
        end
        algp.KRRcost(itr) = sum(abs(saverflast(:) - rf(:)).^2)./2 + ((algp.bSARperc)./2).*sum(abs(APOCS(:)).^2);
        
        % update target phase patterns and calculate error
        %         rfBest = median(sqz(rf(:,raIdx,:)),2);
        %         meanerr=mean(errAll(raIdx,:),2);
        for slIdx = 1:Nsl
            m = A{slIdx}*rf(:,raIdx,slIdx);
            %             if holeReg && any(abs(m) < holeThresh) && meanerr < 2.5
            %                 m = A{slIdx}*rfBest;
            %                 algp.holesExist = true;
            %             elseif any(abs(m) < holeThresh)
            %                 algp.holesExist = true;
            %             end
            dphs{slIdx} = angle(m);
            % save error
            errAll(raIdx,slIdx) = 1/2*norm(m-exp(1i*dphs{slIdx}))^2;
            %             powAll(raIdx,slIdx) = 1/2*beta*norm(rf(:,raIdx,slIdx))^2;
            vopAll(raIdx,slIdx) = 1/2*algp.betaSAR(slIdx)*norm((maps.C{slIdx}*(size(maps.C{slIdx},1)/Nc))*rf(:,raIdx,slIdx))^2;
            % save shimmed pattern
            mAll(:,:,raIdx,slIdx) = embed(m,logical(maps.mask(:,:,slIdx)));
        end
        itr = itr + 1;
        if rem(itr-1,10)==0
            cost(end+1)=mean(errAll(raIdx,:),2);
        end
        if algp.kk==1 && rem(itr-1,10)==0
            %for the first crossval set, keep track of learning curve.
            valcost(end+1)=calcVALcost(rf,Nc,algp);
            disp(valcost(end));
            if algp.cutoff && ((cost(end)>cost(end-1)) || (valcost(end)>valcost(end-1)))
                flag=0; %cut off if either cost increases
            end
            %also exit if the change in cost is less than a threshold.
            diffcost=diff(cost);
            if abs(diffcost(diffcost<0))<=1e-6
                flag=0;
            end
            mindiff=min(abs(diffcost(diffcost<0)));
%         elseif algp.kk~=1 %quit based on # of iters that we needed for first set.
%             if ~exist('iterstorun') || iterstorun.itr==Inf
%                 if exist(['valcost' algp.fstringItr ],'file')
%                     iterstorun=load(['valcost' algp.fstringItr],'itr');
%                 else
%                     iterstorun.itr=Inf;
%                 end
%             end
%             if iterstorun.itr<=itr
%                 flag=0;
%             end
        elseif algp.kk~=1 %quit based on final DIFF of cost that we needed for first set.
            %i.e. if cost is changing by smaller amt than 1st k-fold was at
            %the end, cut it off.
            if ~exist('difftorun') || difftorun.mindiff==0
                if exist(['valcost' algp.fstringItr ],'file')
                    pause(120); %wait to make sure this thing is fully saved first.
                    difftorun=load(['valcost' algp.fstringItr],'mindiff');
                    if ~isfield(difftorun,'mindiff')
                        difftorun.mindiff=0;
                    end
                else
                    difftorun.mindiff=0;
                end
            end
            

            % turn this off since we monitor validation cost.
%             if algp.cutoff && length(cost)>=2 && ((cost(end)>cost(end-1))) %or if cost increases.
%                 flag=0;
%             end
            diffcost=diff(cost);
            mindiff=min(abs(diffcost(diffcost<0)));
            if difftorun.mindiff>=mindiff %quit when cost isn't changing significantly.
                flag=0;
            end
        end
        
        % check stopping criterion
        if algp.cutoff && (norm(sqz(rf(:,raIdx,:))-rfOld)/norm(rfOld) < tol && ~algp.holesExist)
            flag =0;
        elseif itr > algp.maxIters
            flag = 0;
        end
        
        fprintf('Initialization %d. Iteration %d. Err all %f. SAR all %f.\n',raIdx,itr,sum(errAll(raIdx,:),2),sum(vopAll(raIdx,:),2));
        if algp.dofigs
            figure(1);clf;im(mAll(:,:,raIdx,:));
            drawnow
        end
        
        if ~flag
            save(['valcost' algp.fstring],'valcost','itr','cost','mindiff')
        end
        toc
    end % while loop over all updates
    
    fprintf('Final stats for Initialization %d: Iterations %d. Err all %f. SAR all %f.\n',raIdx,itr,sum(errAll(raIdx,:),2),sum(vopAll(raIdx,:),2));
    
end % loop over random starts


function cost=calcVALcost(rf,Nc,algp)
%% calculate validation cost
% get prediction matrix back and apply to get the test RF
rfPOCS = sqz(rf);
% get A matrix
APOCS = [];
testInds=algp.testInds;
for ii = 1:Nc
    %         APOCS = [APOCS; (algp.Finv*real(rfPOCS(ii,:)).').'; (algp.Finv*imag(rfPOCS(ii,:)).').'];
    APOCS(ii,:) = (algp.Finv*rfPOCS(ii,:).');
end
cost=[];
for ii = 1:size(algp.testFeaturesPOCSshim,2)
    Nslpersub=31;
    predPOCS = APOCS*algp.testFeaturesPOCSshim(:,ii);
    %         predPOCS = predPOCS(1:2:2*Nc) + 1i*predPOCS(2:2:2*Nc);
    % get subject and slice indices
    subjInd = floor(testInds(ii)/Nslpersub)+1;
    slInd = rem(testInds(ii),Nslpersub);
    if slInd == 0;slInd = Nslpersub;subjInd = subjInd - 1;end;
    
    tmp = zeros(size(algp.maps{subjInd}.b1,1),size(algp.maps{subjInd}.b1,2));
    for jj = 1:Nc
        tmp = tmp + algp.maps{subjInd}.b1(:,:,slInd,jj)*predPOCS(jj);
    end
    %         cost(ii) = 100*std(abs(tmp(algp.maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(algp.maps{subjInd}.mask(:,:,slInd))));
    cost(ii) = (1/2)*norm(tmp(algp.maps{subjInd}.mask(:,:,slInd))/mean(abs(tmp(algp.maps{subjInd}.mask(:,:,slInd))))-exp(1i.*angle(tmp(algp.maps{subjInd}.mask(:,:,slInd)))))^2;
end
cost=mean(cost);
clearvars -EXCEPT cost



