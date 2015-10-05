function y = MeasurementMatrix_t2s(alpha,f,L,m,n,Psi,mode)
if ~isa(Psi, 'function_handle')
    Psi = @(x) Psi*x;
end


if ~isa(Psi, 'function_handle')
    Psi = @(x) Psi*x;
end
alpha = reshape(alpha,m,n);
if strcmp(mode,'DCT')
    x = idct2(alpha);
else
    x = midwt(alpha,f,log2(m));
end
x = x(:);

y = Psi(x);


