function [x_denoised, r_tilde_num] = SURE_denoise(x_rec_sub, xi_sub, n_sub)

delta_N = 0.05; kappa = 0.05; alpha = 0.1; beta = 0.3; ell_zero = 20; Flag = 1;
ITER = 30;
r_tilde = @(tau) calcu_r_tilde(tau,x_rec_sub , xi_sub, n_sub);
mu_tilde = abs(1/n_sub*norm(x_rec_sub)^2 - xi_sub);
[tau_new] = appro_gd(10*delta_N, kappa, alpha, beta, mu_tilde, ITER, ell_zero, Flag, n_sub, r_tilde);
x_denoised = wthresh(x_rec_sub,'s',tau_new);
r_tilde_num = abs(r_tilde(tau_new));