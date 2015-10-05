%Generates state evolution and compares to observed MSEs for MR-AMP algorithms
clear all;
addpath(genpath('..'));

%Parameters
Mode='Higher';
SamplingRate=.1;
noise_sig=0;
imsize=128;
filename='barbara.png';
iters=30;
N=5;%Number of tests to run.  Must be at least 2
factor =2; % downsampling factor
Im=double(imread(filename));
x_0=imresize(Im,imsize/size(Im,1));
[height, width]=size(x_0);
n=length(x_0(:));
m=round(n*SamplingRate);
Up_matrix = zeros(imsize,imsize/factor);
%downFactor = N/n;
for i = 1:size(Up_matrix,2)
    Up_matrix((i-1)*factor+1:i*factor,i) = 1;
end
%Generate State Evolution and compute intermediate MSEs of D-AMP and D-IT
scale = 2;
True_AMP_MSE_array=zeros(N,iters);
True_DIT_MSE_array=zeros(N,iters);
Predicted_MSE_array=zeros(N,iters);
True_AMP_MSE_array_LR = zeros(N,iters);
Predicted_MSE_array_LR = zeros(N,iters);
x_LR_0 = imresize(x_0,1/factor,'bicubic'); % the Reference LR image
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
    x_LR_t = zeros(n/factor^2,1);
    z_t=y;
    z_LR_t = y;
    x_HR = imresize(x_LR_0,factor,'bicubic');
    noise_sig_LR = sqrt(mean((y-Psi(x_HR(:))).^2));
    Mode = 'Higher';
    A_LR = @(alpha)MeasurementMatrix_t2s_LR_v2(alpha,imsize/factor,imsize/factor,imsize,imsize,Psi,Mode,Up_matrix,scale);
    AT_LR = @(y)MeasurementMatrixTrans_s2t_LR_v2(y,imsize/factor,imsize/factor,imsize,imsize,Psi_T,Mode,Up_matrix,scale);
  

    for iter=2:iters
        factor = 1; % for HR recovery
        scale = 1; % for HR recovery
        [x_tplus1,z_tplus1,pseudo_data] = TVAMP_oneIter(y,x_t,z_t,Psi,Psi_T,imsize,imsize,n,factor,Mode,scale);
        z_t=z_tplus1; x_t=x_tplus1;
        True_AMP_MSE_array(i,iter)=mean((x_0(:)-x_tplus1(:)).^2);
       [ Predicted_MSE_array(i,iter) ] = TVAMP_SE_Prediction( x_0(:), Predicted_MSE_array(i,iter-1), m,n,noise_sig, factor,Mode,scale);
       
        factor = 2; % for HR recovery
        scale = 2.68; % for HR recovery
        [x_LR_tplus1,z_LR_tplus1,pseudo_data] = TVAMP_oneIter(y,x_LR_t,z_LR_t,A_LR,AT_LR,imsize/factor,imsize/factor,n/factor^2,factor,Mode,scale);
        z_LR_t=z_LR_tplus1; x_LR_t=x_LR_tplus1;
        True_AMP_MSE_array_LR(i,iter)=mean((x_LR_0(:)-x_LR_tplus1(:)/scale).^2);
        mse_tmp = TVAMP_SE_Prediction( x_LR_0(:)*scale, Predicted_MSE_array_LR(i,iter-1)*scale^2, m,n/factor^2,noise_sig_LR, factor,Mode,scale);
        Predicted_MSE_array_LR(i,iter)=mse_tmp/scale^2;
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
xlabel('Iteration');
ylabel('MSE');

