function [ MSE_tplus1 ] = TVAMP_SE_Prediction( x_0, MSE_t, m,n,noise_sig, factor,Mode,scale)
%[ MSE_tplus1 ] = DAMP_SE_Prediction( x_0, MSE_t, m,n,noise_sig,denoiser,width,height)
%   This function computes the state evolution prediction for a given
%   denoiser, noise level, and sampling rate
% Input:
%       x_0       : The true signal
%       MSE_t     : The current state evolution MSE esimate
%       m         : number of measurements
%       n         : the signal length
%       noise_sig : the standard deviation of the measurement noise
%       width     : width of the sampled signal
%       height    : height of the sampeled signal. height=1 for 1D signals
%       denoiser  : string that determines which denosier to use. e.g.
%       denoiser='BM3D'
%Output:
%       MSE_tplus1: The next state evolution MSE esimate
%    denoi=@(noisy,sigma_hat) denoise(noisy,sigma_hat,width,height,denoiser);

% default configuration of 'auto' mode for parameterless AMP
if strcmp(Mode,'Single')
    factor = factor;
else
    factor = scale;
end

d=m/n;
Z=randn(n,1);
netsig=sqrt((1/d)*(MSE_t)+ noise_sig^2);
noisy=x_0+netsig*Z;
noisy_norm = noisy/factor;
xi_norm = netsig/factor;
% xi_norm = sqrt(mse_true(i));
x_t_norm = tvdenoise_adapt(noisy_norm,xi_norm);
%  x_t_norm = tvmm_a(pseudo_data_norm,1,(xi_norm));
denoised = x_t_norm*factor;
MSE_tplus1 = mean((denoised(:)-x_0(:)).^2);
end

