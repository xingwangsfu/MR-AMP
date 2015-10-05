function y = MeasurementMatrix_t2s_LR_hybrid(alpha,f,L,m,n,M,N,Psi,mode,alpha_true)


if strcmp(mode,'DCT')
    factor = M/m;
    alpha = reshape(alpha,m,n);
    x_LR = idct2(alpha)/factor;
    It=twotime_1dfilter(x_LR,1,1);
    It=unified_iterative(x_LR, It, 8, 'overlap_fit', 1, 'no_fit', 1);
    %     filter_coef = [1 -5 20 20 -5 1]/32;
    %     x_HR_Wiener = my_interp(x_LR,filter_coef);
    alpha_HR = dct2(It);
    alpha_HR(1:m,1:n) = alpha;
    x = idct2(alpha_HR);
else
    f0 = f{1};
    f1 = f{2};
    alpha = reshape(alpha,m,n);
    alpha_HR = zeros(M,N);
    alpha_HR(1:m,1:n) = alpha;
    x = idwt2d(alpha_HR,f0,f1,L);
end
x = x(:);

y = Psi*x;