% this file is to do multi-resolution compressed sensing
% in spatial domain with total-variation denoiser

clear all;
addpath('C:\Users\xingw\Desktop\research file\xingw\scalable CS');
% input image
img = imread('lena.png');
% img = rgb2gray(img);
img = double(img);
n=128;
[height,width] = size(img);
X = imresize(img,n/height); % the true HR image, resized to 128*128
N = n*n;
factor =2; % downsampling factor
n_LR = floor(n/factor);
Up_matrix = zeros(n,n_LR); % upsampling matrix (prolongation operator)
if mod(factor,1) ==1
    for i = 1:size(Up_matrix,2)
        Up_matrix((i-1)*factor+1:i*factor,i) = 1;
    end
end


delta_cand = [0.2]; % undersampling pool
X_vec = X(:); % vectorize the input image
for i = 1:length(delta_cand)
    
    delta = delta_cand(i);
    M = floor(N*delta);
    
    for index = 1:20
        % Measurement matrix Psi
        Phi = randn(M,N);
        normalized_vector = sqrt(sum(Phi.*Phi,1));
        normalized_matrix = repmat(normalized_vector,M,1);
        Phi = Phi./normalized_matrix;
        
        if ~isa(Phi, 'function_handle')
            Phi_T = @(x) Phi'*x;
            Phi = @(x) Phi*x;
        end
        Y = Phi(X_vec);

        Mode = 'TV'; % For HR image reconstruction, no need to input Mode and scale.
        scale = 1;
        
        % recover the HR image, the output is recovered HR image
        [alpha_rec, nmse, mse_true, mse_pred] = TVAMP(Y, Phi, Phi_T, n, n, N,X,1,Mode,scale);
        alpha_H2L = imresize(alpha_rec,1/factor,'bicubic'); % HRL
        x_LR = imresize(X,1/factor,'bicubic'); % reference LR image with bicubicu interpolator
        x_LR_single = X(1:factor:end,1:factor:end); % reference LR image with prolongation operator
        PSNR_H2L(i,index) =  20*log10(255/sqrt(norm(alpha_H2L(:)-x_LR(:))^2/(N/factor^2)));
        PSNR_HR(i,index) = 20*log10(255/sqrt(norm(alpha_rec(:)-X(:))^2/(N)));
        
        % recover the LR image, the output is the recovered LR image
        Mode = 'Higher'; % two modes, "Single" order for prolongation, and "Higher" for bicubic
        scale = 2.68; % scaling factor for LR image in "Higher" mode, scale=2.68 for factor=2, scale=5 for factor=4
                      %  scale = factor for LR image in "Single" mode
        
        % Construct corresponding new Measurement matrix for LR
        A_LR = @(alpha)MeasurementMatrix_t2s_LR_v2(alpha,n_LR,n_LR,n,n,Phi,Mode,Up_matrix,scale);
        AT_LR = @(y)MeasurementMatrixTrans_s2t_LR_v2(y,n_LR,n_LR,n,n,Phi_T,Mode,Up_matrix,scale);
        [alpha_rec_LR, nmse_LR_tmp, mse_true_LR, mse_pred_LR] = TVAMP(Y, A_LR, AT_LR, n_LR, n_LR, n_LR^2,x_LR,factor,Mode,scale);
        if strcmp(Mode,'Single')
            PSNR_LR(i,index) =  20*log10(255/sqrt(norm(alpha_rec_LR(:)/(factor)-x_LR(:))^2/(N/factor^2)));
        else
            PSNR_LR(i,index) =  20*log10(255/sqrt(norm(alpha_rec_LR(:)/(scale)-x_LR(:))^2/(N/factor^2)));
        end
    end
end

