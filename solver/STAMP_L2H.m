function [x_L2H] = STAMP_L2H(y,A,AT,n,lambda,x_SI,sigma_SI)

% input: y: measurements
%        A: measurement matrix
%        AT: transpose of measurement matrix
%        n: dimension of x
%        lambda: the thresholding parameter in soft-thresholding
%        Params:   Params.x_SI, side information
%                  Params.T, number of iterations
%                  Params.tol, tolerance
%                  Params.sigma_SI: variance of side information
%                  Params.mode: two modes here: auto and minimax, auto for parameterless AMP, minimax for DMM
%                  Params.SIorNot : two modes: SI exists or not
% Output: x_t: reconstructed signal
%         r_tilde_real: the estimated MSE predicted by SURE (Stein's
%         unbiased risk estimate) used in parameterless AMP

% m, n, x_SI, sigma, sigma_SI, alpha_num, epsilon, T, tol, delta, x, mode, SIorNot
% m = Params.m;
% n = Params.n;

% test the number of input
% [m, n] = size(A);
% initial estimate

if ~isa(A, 'function_handle')
    AT = @(x) A'*x;
    A = @(x) A*x;
end

x_t = zeros(n,1);
blksize_H = sqrt(n);
n_L = numel(x_SI);
blksize_L = sqrt(n_L);
r_t = y-A(x_t);
m=length(y);

% default configuration of 'auto' mode for parameterless AMP
T =30;
r_tilde_real = zeros(T,1);


for i = 1:T
    res_energy = norm(r_t)^2/m;
    mu_L = res_energy./sigma_SI;
    xi_L = sigma_SI*res_energy/(sigma_SI+res_energy);
    xi = res_energy;
    pseudo_data = x_t + AT(r_t);
    pseudo_data_mtx = reshape(pseudo_data,blksize_H,blksize_H);
    pseudo_data_L = pseudo_data_mtx(1:blksize_L,1:blksize_L);
    pseudo_data_H_Fhalf = pseudo_data_mtx(blksize_L+1:end,1:blksize_L);
    pseudo_data_H_Shalf = pseudo_data_mtx(:,blksize_L+1:end);
    pseudo_data_H = [pseudo_data_H_Fhalf(:);pseudo_data_H_Shalf(:)];
    pseudo_data_L = mu_L*x_SI(:)/(1+mu_L) + (pseudo_data_L(:))/(1+mu_L);
    [x_t_L, r_tilde_num_L] = SURE_denoise(pseudo_data_L(:), xi_L, n_L);
    [x_t_H, r_tilde_num_H] = SURE_denoise(pseudo_data_H(:), xi, n-n_L);
    x_t = zeros(blksize_H,blksize_H);
    x_t(1:blksize_L,1:blksize_L) = reshape(x_t_L,blksize_L,blksize_L);
    x_t(blksize_L+1:end,1:blksize_L) = reshape(x_t_H(1:(blksize_H-blksize_L)*blksize_L), blksize_H-blksize_L, blksize_L);
    x_t(:,blksize_L+1:end) = reshape(x_t_H((blksize_H-blksize_L)*blksize_L+1:end),blksize_H,blksize_H-blksize_L);
    x_t = x_t(:); 
    r_t = y-A(x_t)+r_t*1/m*(sum(x_t_L~=0)/(1+mu_L)+sum(x_t_H~=0));
end

x_L2H = x_t;






