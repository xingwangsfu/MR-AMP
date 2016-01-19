function y = Phi_LR(x,m,n,M,N,Psi,mode,Up_matrix,scale)

if ~isa(Psi, 'function_handle')
    Psi = @(x) Psi*x;
   % A = @(x) A*x;
end

x = reshape(x,m,n);
%Up_matrix = zeros(N,n);
downFactor = N/n;

if strcmp(mode,'Single')
    x_HR = 1/downFactor*Up_matrix*x*Up_matrix';
else
    %  if downFactor ==2
    x_HR = imresize(x,downFactor,'bicubic')/(scale);
end
y = Psi(x_HR(:));