function [x_tplus1,z_tplus1,pseudo_data] = AMP_oneIter(y, x_t, z_t, A, AT, n, Mode)

if ~isa(A, 'function_handle')
    AT = @(x) A'*x;
    A = @(x) A*x;
end

if nargin < 6
    error('Wrong Number of Input Parameters!');
else if nargin < 7
        Mode = 'Minimax';
    end
end
m=length(y);
delta = m/n;
if strcmp(Mode,'Auto')
    % default configuration of 'auto' mode for parameterless AMP
    res_energ = norm(z_t)^2/m;
    xi = res_energ;
    pseudo_data = x_t + AT(z_t);
    [x_tplus1, r_tilde_num] = SURE_denoise(pseudo_data, xi, n);
    z_tplus1 = y-A(x_tplus1)+z_t*sum(x_tplus1~=0)/(m);
else
    [alpha,lambda] = optimal_threshold_DMM(delta);
    res_energ = norm(z_t)^2/m;
    xi = res_energ;
    pseudo_data = x_t + AT(z_t);
    x_tplus1 = soft_threshold(pseudo_data, lambda*sqrt(xi));
    z_tplus1 = y-A(x_tplus1)+z_t*sum(x_tplus1~=0)/(m);
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
