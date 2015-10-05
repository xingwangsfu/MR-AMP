function [ MSE_tplus1 ] = AMP_SE_Prediction( x_0, MSE_t, m,n,noise_sig, Mode)
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

if nargin < 5
    error('Wrong Number of Input Parameters!');
else if nargin < 6
        Mode = 'Minimax';
    end
end
d=m/n;
Z=randn(n,1);
netsig=sqrt((1/d)*(MSE_t)+ noise_sig^2);
noisy=x_0+netsig*Z;

if strcmp(Mode,'Auto')
    % default configuration of 'auto' mode for parameterless AMP
    [denoised, r_tilde_num] = SURE_denoise(noisy, netsig^2, n);
    MSE_tplus1=sum((denoised-x_0).^2)/n;
else
    [alpha,lambda] = optimal_threshold_DMM(d);
    denoised = soft_threshold(noisy, lambda*netsig);
    MSE_tplus1=sum((denoised-x_0).^2)/n;
end
end

function x_thresh = soft_threshold(x_unthresh, theta)
        
    if isscalar(theta)
            x_thresh = wthresh(x_unthresh, 's', theta);
    else
            x_thresh = zeros(size(x_unthresh));
            index_1 = find(x_unthresh>theta);
            index_2 = find(abs(x_unthresh)<=theta);
            index_3 = find(x_unthresh<(-1*theta));
            x_thresh(index_1) = x_unthresh(index_1) - theta(index_1);
            x_thresh(index_2) = 0;
            x_thresh(index_3) = x_unthresh(index_3)+theta(index_3);
    end
end
