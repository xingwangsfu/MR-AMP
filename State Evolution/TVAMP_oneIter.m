function [x_tplus1,z_tplus1,pseudo_data] = TVAMP_oneIter(y,x_t,z_t,A,AT,m,n,N,factor,Mode,scale)

if ~isa(A, 'function_handle')
    AT = @(x)A'*x;
    A = @(x)A*x;
end

% x_t = zeros(N,1);
% z_t = y;
if strcmp(Mode,'Single')
    factor = factor;
else
    factor = scale;
end

M=length(y);
delta = M/N;
% [alpha,lambda] = optimal_threshold_TV(delta); % the remaining problem is to determine lambda given delta
T = 30;
%lambda = 20;
tol = 1e-6;
%opts.lambda = lambda;
%B = speye(N,N);
pseudo_data = AT(z_t)+x_t;
xi = sqrt(norm(z_t)^2/M);
% mse_true(i) = norm(pseudo_data(:)/factor-X(:))^2/N;
% mse_pred(i) = xi^2/factor^2;

% normalize
pseudo_data_norm = pseudo_data/factor;
xi_norm = xi/factor;
% xi_norm = sqrt(mse_true(i));
x_t_norm = tvdenoise_adapt(pseudo_data_norm,xi_norm);
%  x_t_norm = tvmm_a(pseudo_data_norm,1,(xi_norm));
x_tplus1 = x_t_norm*factor;
% nmse(i) =  norm(x_t(:)/factor-X(:))^2/norm(X(:))^2;
eta=randn(1,N);
epsilon=(max(pseudo_data(:))/1000+eps);
x_rec_MC_norm = tvdenoise_adapt(pseudo_data_norm+epsilon/factor*eta',xi_norm);
%    x_rec_MC_norm = tvmm_a(pseudo_data_norm+epsilon/factor*eta',1,(xi_norm));
x_rec_MC = x_rec_MC_norm*factor;
div=eta*((x_rec_MC-x_tplus1)/epsilon); % numerical monte carlo to compute the gradient
z_tplus1=y-A(x_tplus1)+1/M.*z_t.*div;
%     if (norm(y-A(x_t))/norm(y)<tol)
%         break;
%     end


