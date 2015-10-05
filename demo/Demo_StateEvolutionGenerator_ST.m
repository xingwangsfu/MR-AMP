%Generates state evolution and compares to observed MSEs for MR-AMP algorithms
clear all;
addpath(genpath('..'));

%Parameters
Mode='DCT';
SamplingRate=.2;
noise_sig=0;
imsize=128;
filename='barbara.png';
iters=30;
N=5;%Number of tests to run.  Must be at least 2
factor =2;
Im=double(imread(filename));
x_0=imresize(Im,imsize/size(Im,1));
[height, width]=size(x_0);
n=length(x_0(:));
m=round(n*SamplingRate);
%Generate State Evolution and compute intermediate MSEs of D-AMP and D-IT
if strcmp(Mode,'DCT')
    alpha_0= dct2(x_0);
    h = [];
    L =[];
    f = [];
else
    % 9-7 filter
    [h0,h1,f0,f1] = filter9_7();
    L = floor(log2(imsize))-3;
    alpha_0 = dwt2d(x_0,h0,h1,L);
    h{1} = h0;
    h{2} = h1;
    f{1} = f0;
    f{2} = f1;
end
alpha_LR_0 = alpha_0(1:imsize/factor,1:imsize/factor); 
True_AMP_MSE_array=zeros(N,iters);
True_DIT_MSE_array=zeros(N,iters);
Predicted_MSE_array=zeros(N,iters);
True_AMP_MSE_array_LR = zeros(N,iters);
Predicted_MSE_array_LR = zeros(N,iters);
% x_LR_0 = imresize(x_0,1/factor,'bicubic');
x_LR_0 = 1/factor*idct2(alpha_LR_0);
True_AMP_MSE_array(:,1)=mean(x_0(:).^2);
Predicted_MSE_array(:,1)=mean(x_0(:).^2);
True_AMP_MSE_array_LR(:,1) = mean(x_LR_0(:).^2);
Predicted_MSE_array_LR(:,1) = mean(x_LR_0(:).^2);
for i=1:N
    noise_sig = 0;
    M=randn(m,n);
    for j = 1:n
        M(:,j) = M(:,j) ./ sqrt(sum(abs(M(:,j)).^2));
    end
    
    if ~isa(M, 'function_handle')
        Psi_T = @(x) M'*x;
        Psi = @(x) M*x;
    end
    
    y=Psi(x_0(:))+noise_sig*randn(m,1);
    x_t=zeros(n,1);
    x_LR_t=zeros(n/factor^2,1);
    z_t=y;
    z_LR_t=y;
    alpha_HR = zeros(imsize,imsize);
    alpha_HR(1:imsize/factor,1:imsize/factor) = alpha_LR_0;
    x_HR = idct2(alpha_HR);
    noise_sig_LR = sqrt(mean((y-Psi(x_HR(:))).^2)); % the new noise level: Gaussian noise and LR approximation error
    % DCT or Wavelet
    A = @(alpha)MeasurementMatrix_t2s(alpha,f,L,imsize,imsize,Psi,Mode);
    AT = @(y)MeasurementMatrixTrans_s2t(y,h,L,imsize,imsize,Psi_T,Mode);
    A_LR = @(alpha)MeasurementMatrix_t2s_LR(alpha,f,L,imsize/factor,imsize/factor,imsize,imsize,Psi,Mode);
    AT_LR = @(y)MeasurementMatrixTrans_s2t_LR(y,h,L,imsize/factor,imsize/factor,imsize,imsize,Psi_T,Mode);
    
    for iter=2:iters
        [x_tplus1,z_tplus1,pseudo_data] = AMP_oneIter(y, x_t, z_t, A, AT, n);
        z_t=z_tplus1; x_t=x_tplus1;
        True_AMP_MSE_array(i,iter)=mean((alpha_0(:)-x_tplus1(:)).^2);
        Predicted_MSE_array(i,iter) = AMP_SE_Prediction(alpha_0(:), Predicted_MSE_array(i,iter-1), m,n,noise_sig);  
        
        [x_LR_tplus1,z_LR_tplus1,pseudo_data] = AMP_oneIter(y, x_LR_t, z_LR_t, A_LR, AT_LR, n/factor^2);
        z_LR_t=z_LR_tplus1; x_LR_t=x_LR_tplus1;
        True_AMP_MSE_array_LR(i,iter)=mean((alpha_LR_0(:)/factor-x_LR_tplus1(:)/factor).^2);
        mse_tmp = AMP_SE_Prediction(alpha_LR_0(:), Predicted_MSE_array_LR(i,iter-1)*factor^2, m,n/factor^2,noise_sig_LR);
        Predicted_MSE_array_LR(i,iter)=mse_tmp/factor^2;
    end
end

% MSE of AMP
True_AMP_MSE=mean(True_AMP_MSE_array)';
True_AMP_MSE_std=std(True_AMP_MSE_array)';
Predicted_MSE=mean(Predicted_MSE_array)';
Predicted_MSE_std=std(Predicted_MSE_array)';

% MSE of LR-AMP
True_AMP_MSE_LR=mean(True_AMP_MSE_array_LR)';
True_AMP_MSE_std_LR=std(True_AMP_MSE_array_LR)';
Predicted_MSE_LR=mean(Predicted_MSE_array_LR)';
Predicted_MSE_std_LR=std(Predicted_MSE_array_LR)';
%Plot Results
h=figure;
hold
errorbar(0:29,Predicted_MSE,Predicted_MSE_std,'-.b');
errorbar(0:29,True_AMP_MSE,True_AMP_MSE_std,'--g');
xlabel('Iteration');
ylabel('MSE');

h=figure;
hold
errorbar(0:29,Predicted_MSE_LR,Predicted_MSE_std_LR,'-.b');
errorbar(0:29,True_AMP_MSE_LR,True_AMP_MSE_std_LR,'--g');
title([denoiser,'-AMP and ']);
xlabel('Iteration');
ylabel('MSE');


