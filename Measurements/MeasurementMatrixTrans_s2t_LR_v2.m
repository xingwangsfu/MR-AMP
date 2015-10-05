function x = MeasurementMatrixTrans_s2t_LR_v2(y,m,n,M,N,Psi_T,mode,Up_matrix, scale)


if ~isa(Psi_T, 'function_handle')
    Psi_T = @(x) Psi_T*x;
end
tmp = Psi_T(y);

tmp_2d = reshape(tmp,M,N);
downFactor = N/n;
if strcmp(mode,'Single')
    
    x = 1/downFactor*Up_matrix'*tmp_2d*Up_matrix;
else
    %  if downFactor ==2
    x = imresize(tmp_2d,1/downFactor,'bicubic')*(scale);
end
x = x(:);
