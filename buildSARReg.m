function C = buildSARReg(vops,Np)

% builds an SAR regularization matrix using a set of VOPs,
% replicated for Np pulses. The output matrix has column dimension
% organized coil major, time minor.

[Nc,~,Nvops] = size(vops);

% factor the vops
Cind = zeros(size(vops));
for ii = 1:Nvops
    [v,d] = eig(vops(:,:,ii));
    d(d < 0) = 0;
    Cind(:,:,ii) = sqrt(d)*v';
end

% reshape and duplicate to Np pulses
C = sparse(Nc*Nvops*Np,Nc*Np);
for ii = 1:Nc
    Cvec = squeeze(Cind(:,ii,:)); % extract all VOP columns
    Cvec = Cvec(:); % columnize
    % duplicate to Npulses, build block diagonal matrix
    % and concatenate into C matrix with the rest of the Tx channels
    C(:,(ii-1)*Np+1:ii*Np) = kron(speye(Np),Cvec);
end
