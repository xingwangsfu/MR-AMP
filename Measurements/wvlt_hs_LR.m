% wvlt_hs.m
%
% 3D wavelet transform for hyperspectral imaging
%
% Usage: w = wvlt_hs(x, h)
%
% x - input signal of size YxMxN, where MxN is the image size and Y is the
% number of spectral bands
%
% h - Wavelet filter name
%
% w - 3D wavelet coefficient vector for x
%
% Written by: Marco F. Duarte, Rice University
% Created: November 13 2006

function [w_LR] = wvlt_hs_LR(x, h, factor)

factor_t = factor(1); % downsampling factor in the temporal domain
factor_r = factor(2); % downsampling factor in the spatial domain

N = size(x,2);
Y = size(x,1);

Y_LR = Y/factor_t;
N_LR = N/factor_r;
if Y > 1,
    wvlts = zeros(size(x));
    for i=1:N,
        for j=1:N,
            % Perform 1D DWT (spectra) on pixel (i,j)
            wvlts(:,i,j) = mdwt(x(:,i,j),h,log2(Y));
        end
    end
    wvlts_LR = wvlts(1:Y_LR,:,:);
else
    wvlts = x;
end

w_LR = zeros([Y_LR N_LR N_LR]);
for i=1:Y_LR,
    % Perform 2D DWT on snapshot i (space)
    tmp = reshape(mdwt(squeeze(wvlts_LR(i,:,:)),h,log2(N)),[1 N N]);
    w_LR(i,:,:) = tmp(1,1:N_LR,1:N_LR);
end
w_LR = w_LR(:);
%w = w(:);

