% This function fourierShapeDescriptor calculates fourier shape descriptors
% of simple binary region boundary according to the method described by Rafael C. 
% Gonzalez, Richard E. Woods "Digital Image Processing" 3rd ed. Chapter 11 (p840)
%
% B = fourierShapeDescriptor(A,r)
%
% inputs,
%   A : The binary image containing the region of interest. This program
%       was tested for single connected component (CC) only. Multiple CCs will
%       cause only one CC to be described.
%   r : Ratio of boundary pixel to fourier coefficients. Large ratio means
%       the shape would be described by less fourier coeffiecient. r = 1
%       means will yield the unapproximated region boundary in A.
%       Typically, r could be less than 100 without sacrificing boundary
%       quality.
% ***JDI edited r to be the exact number of fourier descriptors requested.
% outputs,
%   B: Output binary image reconstructed by the number of fourier
%   coefficients specified by the user.
%
%
% Written by Weicheng Kuo, 07/10/2014

function [B,fdes]=fourierShapeDescriptor(A,r)

bdry_idx = bwboundaries(A);
%***JDI: find max length (outer) boundary:
for ii=1:size(bdry_idx,1)
    blength(ii)=size(bdry_idx{ii,1},1);
end
[~,maxind]=max(blength);
idx_x = bdry_idx{maxind,1}(:,1);
idx_y = bdry_idx{maxind,1}(:,2);
idx_xn = (idx_x - mean(idx_x))/std(idx_x);
idx_yn = (idx_y - mean(idx_y))/std(idx_y);
idx = idx_xn + 1i*idx_yn;
N = numel(idx_x);
n = r;%(N/r);
nstart = floor(N/2)-floor(n/2)+1;
nend = floor(N/2)+floor(n/2)+1;
% if nend>N
%     nend = N;
% end
[h,w] = size(A);

fdes = fftshift(fft(idx));
fdes = fdes(nstart:nend);
fdes_lp = ifftshift([zeros(nstart-1,1);fdes;zeros(N-nend,1)]);
fshape_z = ifft(fdes_lp);
fshape_x = real(fshape_z);
fshape_y = imag(fshape_z);
fshape_x = (fshape_x-mean(fshape_x))/std(fshape_x);
fshape_y = (fshape_y-mean(fshape_y))/std(fshape_y);
fshape_x = nearest(mean(idx_x)+fshape_x*std(idx_x));
fshape_y = nearest(mean(idx_y)+fshape_y*std(idx_y));

B = zeros(h,w);
fdes_idx = sub2ind(size(A),fshape_x,fshape_y);
B(fdes_idx) = 1;
