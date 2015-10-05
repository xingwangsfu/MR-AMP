function [x_t, r_tilde_real] = SI_AMP(y, A, AT, n, lambda, Params )
% SI-AMP
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

% test the number of inputs
if nargin < 5
    error('Wrong Number of Input Parameters');
else if nargin < 6 % load the defaul setting of Params
        Params.x_SI = 0; % no side information
        Params.T = 30; % number of iterations
        Params.tol = 1e-8; % tolerance
        Params.sigma_SI = 0;
        Params.mode = 'Auto'; % two modes here: auto and minimax, auto for parameterless AMP, minimax for DMM
        Params.SIorNot = 'Not'; % two modes: SI exists or not
    end
end

x_SI = Params.x_SI;
sigma_SI = Params.sigma_SI;
T = Params.T;
tol = Params.tol;
mode = Params.mode;
SIorNot = Params.SIorNot;
% [m, n] = size(A);
% initial estimate

if ~isa(A, 'function_handle')
    AT = @(x) A'*x;
    A = @(x) A*x;
end

x_t = zeros(n,1);
r_t = y-A(x_t);
m=length(y);

% default configuration of 'auto' mode for parameterless AMP

r_tilde_real = zeros(T,1);

if strcmp(mode, 'Auto') % 'auto' mode, parameterless AMP
    for i = 1:T
        res_energ = norm(r_t)^2/m;
        if ~strcmp(SIorNot, 'Not')
            mu = res_energ./sigma_SI;
            xi = sigma_SI*res_energ/(sigma_SI+res_energ);
        else
            mu = 0;
            xi = res_energ;
        end
        pseudo_data = mu*x_SI/(1+mu) + (x_t + AT(r_t))/(1+mu);
        [x_t, r_tilde_num] = SURE_denoise(pseudo_data, xi, n);
        r_tilde_real(i) = r_tilde_num;
        r_t = y-A(x_t)+r_t*sum(x_t~=0)/(m*(1+mu));
    end
    
else % 'minimax' mode
    for i = 1:T
        res_energ = norm(r_t)^2/m;
        if ~strcmp(SIorNot, 'Not')
            mu = res_energ./sigma_SI;
            xi = sigma_SI*res_energ/(sigma_SI+res_energ);
        else
            mu = 0;
            xi = res_energ;
        end
        pseudo_data = mu*x_SI/(1+mu) + (x_t + AT(r_t))/(1+mu);
        x_t = wthresh(pseudo_data,'s',lambda*sqrt(xi));
        r_t = y-A(x_t)+r_t*sum(x_t~=0)/(m*(1+mu));
    end
end
end






