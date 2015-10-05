function y = blk_mtx_hs_freq(alpha,Phi,h,N,S)

x = iwvlt_hs(alpha,h,S);
y = Phi*x(:);
end
