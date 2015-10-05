% this file is to compare TV-AMP-2D with TVAL3 algorithm

clear all;
addpath('C:\Users\xingw\Desktop\research file\xingw\scalable CS');
% input image
img = imread('barbara.png');
% img = rgb2gray(img);
img = double(img);
[height,width] = size(img);
X = imresize(img,128/height); % the true HR image, resized to 128*128
n=128;
N = n*n;
factor =1.5;
n_LR = floor(n/factor);
Up_matrix = zeros(n,n_LR);
factor = n/n_LR;
%downFactor = N/n;
if mod(factor,1) ==1
    for i = 1:size(Up_matrix,2)
        Up_matrix((i-1)*factor+1:i*factor,i) = 1;
    end
end


delta_cand = [0.2]; % undersampling pool
sigma_20 = 20^2;
sigma_40 = 40^2;

for i = 1:length(delta_cand)
    
    delta = delta_cand(i);
    
    M = floor(N*delta);
    opts_TVAL3.mu =  2^7;
    opts_TVAL3.beta = 2^4;
    opts_TVAL3.mu0 = opts_TVAL3.mu*0.1;
    opts_TVAL3.beta0 = opts_TVAL3.beta*0.1;
    opts_TVAL3.maxit=200;
    opts_TVAL3.maxcnt = 10;
    opts_TVAL3.tol_inn = 1e-3;
    opts_TVAL3.tol = 1e-6;
    opts_TVAL3.init=0;
    % opts_TVAL3.TVnorm = 1; % Anisortopic setting
    X_vec = X(:);
    for index = 1:20
        % Measurement matrix Psi
        Psi = randn(M,N);
        normalized_vector = sqrt(sum(Psi.*Psi,1));
        normalized_matrix = repmat(normalized_vector,M,1);
        Psi = Psi./normalized_matrix;
        
        Y = Psi*X_vec;
        w = sqrt(sigma_20)*randn(M,1);
        Y_noisy_20 = Y+w;
        w = sqrt(sigma_40)*randn(M,1);
        Y_noisy_40 = Y+w;
        
        %         Mode = 'TV';
        %         % Construct corresponding new Measurement matrix for HR and LR
        %         A = @(alpha)MeasurementMatrix_t2s(alpha,h,L,n,n,Psi,Mode);
        %         AT = @(y)MeasurementMatrixTrans_s2t(y,h,L,n,n,Psi,Mode);
        %         scale = 1;
        %
        %         tic;
        %         [alpha_rec, nmse, mse_true, mse_pred] = TVAMP(Y, Psi, Psi', n, n, N,X,1,Mode,scale);
        %         t_HR(i,index) = toc;
        %         alpha_H2L = imresize(alpha_rec,1/factor,'bicubic');
        %         x_LR = imresize(X,1/factor,'bicubic');
        %         PSNR_H2L(i,index) =  20*log10(255/sqrt(norm(alpha_H2L(:)-x_LR(:))^2/(N/factor^2)));
        %         PSNR_HR(i,index) = 20*log10(255/sqrt(norm(alpha_rec(:)-X(:))^2/(N)));
        
        opts.mu = 2^11;
        opts.beta = 2^6;
        opts.mu0 = 2^4;       % trigger continuation shceme
        opts.beta0 = 2^-2;    % trigger continuation shceme
        opts.maxcnt = 10;
        opts.tol_inn = 1e-3;
        opts.tol = 1E-6;
        opts.maxit = 300;
        opts.nonneg = true;
        opts.isreal = true;
        tic;
        [estIm, out] = TVAL3(Psi,Psi',Y,n,n,opts_TVAL3);
        t_HR_TVAL3(i,index) = toc;
        PSNR_HR_TVAL3(i,index) = 20*log10(255/sqrt(norm(estIm(:)-X(:))^2/(N)));
        index
    end
end







% set the optional parameters
