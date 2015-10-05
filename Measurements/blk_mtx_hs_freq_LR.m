function y = blk_mtx_hs_freq_LR(alpha,Phi,h,N,S, factor)

% x = iwvlt_hs(alpha,h,S);
x = iwvlt_hs_LR(alpha, h, N, S, factor);
y = Phi*x(:);
end