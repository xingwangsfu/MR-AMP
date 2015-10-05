% iwvlt_hs.m
%
% 3D inverse wavelet transform for hyperspectral images
%
% Usage: x = iwvlt_hs(w, h, S)
%
% w - input wavelet coefficient vector
%
% h - Wavelet filter
%
% S - number of spectral bands
%
% x - output data cube
%
% Written by: Marco F. Duarte, Rice University
% Created: November 5 2006
%
% Uses Rice Wavelet Toolbox, http://dsp.rice.edu/rwt

function x = iwvlt_hs_LR(w_LR, h, M, S, factor)

factor_t = factor(1);
factor_s = factor(2);
S_LR = S/factor_t;
M_LR = M/factor_s;
w = reshape(w_LR,[S_LR M_LR M_LR]);
wvltx = zeros([S_LR M M]);
for i=1:S_LR,
    % Perform 2D IWT on snapshot i (space)
    tmp_w = zeros([1, M, M]);
    tmp_w(1,1:M_LR,1:M_LR) = w(i,:,:);
    tmp = reshape(midwt(squeeze(tmp_w),h,log2(M)),[1 M M]);
    wvltx(i,:,:) = tmp(1,:,:);
end

if S == 1,
    x = wvltx;
else
    x = zeros([S M M]);
    for i=1:M,
        for j=1:M,
            tmp_w = zeros([S,1,1]);
            tmp_w(1:S_LR,:,:) = wvltx(:,i,j);
            % Perform 1D DWT (spectra) on pixel (i,j)
            x(:,i,j) = midwt(tmp_w,h,log2(S));
        end
    end
end
