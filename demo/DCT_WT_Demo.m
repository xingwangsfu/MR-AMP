% This file is to do multi-resolution compressed sensing
% in transform domain (DCT or wavelet) with soft-thresholding denoiser
clear all;
% input image
img = imread('lena.png');
[height,width,depth] = size(img);
if depth ~=1
    img = rgb2gray(img);
end
img = double(img);
n=128;
X = imresize(img,n/height); % the true HR image, resized to 128*128

Mode = 'DCT'; % two modes are provided: 'DCT', and 'Wavelet'
factor =2; % downsampling factor
if strcmp(Mode,'DCT')% compute the DCT coefficients
    alpha_true = dct2(X);
else %compute the wavelet coefficients
    % 9-7 filter
    h = daubcqf(8);
    f = h;
    L = floor(log2(n));
    alpha_true = mdwt(X,h,L);
end

N = length(alpha_true(:));
n_LR = floor(n/factor);
alpha_true_LR = alpha_true(1:n_LR,1:n_LR); % the DCT (wavelet) coefficient for the LR image
if strcmp(Mode,'DCT')
    h = [];
    f= [];
    L = [];
end

delta_cand = [0.1]; % undersampling rate
X_vec = X(:); % vectorize the input image
for i = 1:length(delta_cand)
    
    delta = delta_cand(i);
    
    M = floor(N*delta);
    
    if (factor == 4) % compute the thresholding parameter of soft-thresholding function using "minimax" rule
        [alpha,lambda] = optimal_threshold_DMM(delta);
        [alpha,lambda_LR] = optimal_threshold_DMM(factor^2*delta);
    else
        lambda = 0;
        lambda_LR = 0;
    end
    
    % Monte Carlo simulation, repeat 20 times
    for index = 1:20
        
        % Generate the measurement matrix
        Phi = randn(M,N);
        normalized_vector = sqrt(sum(Phi.*Phi,1));
        normalized_matrix = repmat(normalized_vector,M,1);
        Phi = Phi./normalized_matrix;
        
        if ~isa(Phi, 'function_handle')
            Phi_T = @(x) Phi'*x;
            Phi = @(x) Phi*x;
        end
        
        Y = Phi(X_vec);
        % Construct corresponding new Measurement matrix for HR and LR
        
        A = @(alpha)Psi_t2s(alpha,f,L,n,n,Phi,Mode);
        AT = @(y)PsiT_s2t(y,h,L,n,n,Phi_T,Mode);
        A_LR = @(alpha)Psi_t2s_LR(alpha,f,L,n/factor,n/factor,n,n,Phi,Mode);
        AT_LR = @(y)PsiT_s2t_LR(y,h,L,n/factor,n/factor,n,n,Phi_T,Mode);
        
        
        % recover the HR image, the output is the DCT coefficients for HR
        % image
        [alpha_rec,xi] = SI_AMP(Y, A, AT, N, lambda); % noiseless case
        
        % recover the LR image, the output is the high-frequency part of
        % DCT coefficients for HR
        [alpha_rec_LR,xi_LR] = SI_AMP(Y, A_LR, AT_LR, n_LR^2, lambda_LR); % noiseless case
        
        %  Compute PSNR and NMSE  for HR and LR
        if strcmp(Mode,'DCT')
            img_rec_LR = 1/factor*idct2(reshape(alpha_rec_LR,n_LR,n_LR)); % the LR image recovered by our proposed method
            img_rec_HR = idct2(reshape(alpha_rec,n,n)); % the recovered HR image
            tmp = reshape(alpha_rec,n,n);
            alpha_H2L = tmp(1:n_LR,1:n_LR);
            img_rec_H2L = 1/factor*idct2(alpha_H2L);
            img_LR = 1/factor*idct2(reshape(alpha_true_LR,n_LR,n_LR)); % the true LR image
        else
            img_rec_LR = 1/factor*midwt(reshape(alpha_rec_LR,n_LR,n_LR),f,L-log2(factor));
            %   img_rec_LR_noisy_40 = 1/factor*idwt2d(reshape(alpha_rec_LR_noisy_40,n/factor,n/factor),f0,f1,L-log2(factor));
            img_rec_HR = midwt(reshape(alpha_rec,n,n),f,L);
            tmp = reshape(alpha_rec,n,n);
            alpha_H2L = tmp(1:n_LR,1:n_LR);
            img_rec_H2L = 1/factor*midwt(alpha_H2L,f,L-log2(factor));
            img_LR = 1/factor*midwt(reshape(alpha_true_LR,n_LR,n_LR),f,L-log2(factor));
        end
        
        PSNR_LR(i,index) = 20*log10(255/sqrt(norm(img_rec_LR(:)-img_LR(:))^2/(prod(size(img_LR)))));
        nmse_LR(i,index) = norm([alpha_rec_LR(:)]-alpha_true_LR(:))^2/norm(alpha_true_LR(:))^2;
        PSNR_HR(i,index) = 20*log10(255/sqrt(norm(img_rec_HR(:)-X(:))^2/(prod(size(X)))));
        nmse_HR(i,index) =  norm(alpha_true(:)-alpha_rec(:))^2/norm(alpha_true(:))^2;
        PSNR_H2L(i,index) = 20*log10(255/sqrt(norm(img_rec_H2L(:)-img_LR(:))^2/(prod(size(img_LR)))));
        nmse_H2L(i,index) = norm([alpha_H2L(:)]-alpha_true_LR(:))^2/norm(alpha_true_LR(:))^2;
    end
end
% clear tmptmp;

