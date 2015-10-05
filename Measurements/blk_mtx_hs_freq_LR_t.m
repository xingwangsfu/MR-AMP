function alpha = blk_mtx_hs_freq_LR_t(y,Phi_T,h,N,S,factor)
x = Phi_T*y;
x = reshape(x,[S N N]);
[alpha] = wvlt_hs_LR(x, h, factor);
% alpha = wvlt_hs(x,h);
end