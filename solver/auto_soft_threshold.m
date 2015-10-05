function [x_denoise] = auto_soft_threshold(x_noisy, sigma)
% default configuration of 'auto' mode for SURE in parameterless AMP
delta_N = 0.05; kappa = 0.05; alpha = 0.1; beta = 0.3; ell_zero = 20; Flag = 1;
ITER = 30;

x_noisy = x_noisy(:);
x_rec_sub = x_noisy;
xi_sub = sigma;
n_sub = prod(size(x_noisy));
r_tilde = @(tau) calcu_r_tilde(tau,x_rec_sub , xi_sub, n_sub);
mu_tilde = abs(1/n_sub*norm(x_rec_sub)^2 - xi_sub);
[tau_new] = appro_gd(10*delta_N, kappa, alpha, beta, mu_tilde, ITER, ell_zero, Flag, n_sub, r_tilde);
x_denoise = wthresh(x_noisy,'s',tau_new);
