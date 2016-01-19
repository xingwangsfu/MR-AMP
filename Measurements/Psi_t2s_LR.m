function y = Psi_t2s_LR(alpha,f,L,m,n,M,N,Psi,mode)

if ~isa(Psi, 'function_handle')
    Psi = @(x) Psi*x;
end


if strcmp(mode,'DCT')
    alpha = reshape(alpha,m,n);
    alpha_HR = zeros(M,N);
    alpha_HR(1:m,1:n) = alpha;
    x = idct2(alpha_HR);
else

    alpha = reshape(alpha,m,n);
    alpha_HR = zeros(M,N);
    alpha_HR(1:m,1:n) = alpha;
    x = midwt(alpha_HR,f,L);
end
x = x(:);

y = Psi(x);