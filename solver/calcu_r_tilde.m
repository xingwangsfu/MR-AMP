function b = calcu_r_tilde(tau, x_rec_sub , xi_sub, n_sub)

x_soft = wthresh(x_rec_sub,'s',tau);

b = (norm(x_soft-x_rec_sub)^2 + n_sub*xi_sub + 2*xi_sub*(sum(x_soft~=0)-n_sub));