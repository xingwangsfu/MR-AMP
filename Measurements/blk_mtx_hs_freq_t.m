function alpha = blk_mtx_hs_freq_t(y,Phi_T,h,N,S)
x = Phi_T*y;
x = reshape(x,[S N N]);
alpha = wvlt_hs(x,h);
end