function alpha_LR = MeasurementMatrixTrans_s2t_LR(y,h,L,m,n,M,N,Psi_T,mode)

if ~isa(Psi_T, 'function_handle')
    Psi_T = @(x) Psi_T*x;
   % A = @(x) A*x;
end

tmp = Psi_T(y);

tmp_2d = reshape(tmp,M,N);
% L = 2;
if strcmp(mode,'DCT')
    alpha_HR = dct2(tmp_2d);
    alpha_LR = alpha_HR(1:m,1:n);
else
    alpha_HR = mdwt(tmp_2d,h,L);
    alpha_LR = alpha_HR(1:m,1:n);
end
alpha_LR = alpha_LR(:);