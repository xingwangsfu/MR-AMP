function [x_t, nmse, x_t_mtx, z_t_mtx] = TVAMP(y,A,AT,m,n,N,X,factor,Mode,scale)

if ~isa(A, 'function_handle')
    AT = @(x)A'*x;
    A = @(x)A*x;
end

x_t = zeros(N,1);
z_t = y;
if strcmp(Mode,'Single')
    factor = factor;
else
    factor = scale;
end

M=length(y);
delta = M/N;
T = 30;
%lambda = 20;
tol = 1e-6;
%opts.lambda = lambda;
%B = speye(N,N);
x_t_mtx = cell(T,1);
z_t_mtx = cell(T,1);
MSE_true(1) = mean((X(:).^2));
for i = 1:T
    pseudo_data = AT(z_t)+x_t;
    xi = sqrt(norm(z_t)^2/M);
    % normalize
    pseudo_data_norm = pseudo_data/factor;
    xi_norm = xi/factor;
    x_t_norm = tvdenoise_adapt(pseudo_data_norm,xi_norm);
    x_t = x_t_norm*factor;
    x_t_mtx{i} = x_t;
    nmse(i) =  norm(x_t(:)/factor-X(:))^2/norm(X(:))^2;
    eta=randn(1,N);
    epsilon=(max(pseudo_data(:))/1000+eps);
    x_rec_MC_norm = tvdenoise_adapt(pseudo_data_norm+epsilon/factor*eta',xi_norm);
    x_rec_MC = x_rec_MC_norm*factor;
    div=eta*((x_rec_MC-x_t)/epsilon); % numerical monte carlo to compute the gradient
    z_t=y-A(x_t)+1/M.*z_t.*div;
    z_t_mtx{i} = z_t;
end
x_t= reshape(x_t,m,n);
