function [alpha] = MeasurementMatrixTrans_s2t(y,h,L,M,N,Psi_T,mode)

if ~isa(Psi_T, 'function_handle')
    Psi_T = @(x) Psi_T*x;
   % A = @(x) A*x;
end

tmp = Psi_T(y);

tmp_2d = reshape(tmp,M,N);
% L=2;

if strcmp(mode,'DCT')
    alpha = dct2(tmp_2d);
else
    alpha = mdwt(tmp_2d,h,log2(N));
  %  [alpha] = dwt2d(tmp_2d,h0,h1,L);
end
alpha = alpha(:);