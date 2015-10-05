function num_change = CountChange(x,m,n)

x = reshape(x,m,n);

[Dx,Dy] = ForwardD(x);

W = sqrt(Dx.*conj(Dx) + Dy.*conj(Dy));

W = W(:);

num_change = length(W~=0);





function [Dux,Duy] = ForwardD(U)
% [ux,uy] = D u

Dux = [diff(U,1,2), U(:,1) - U(:,end)];
Duy = [diff(U,1,1); U(1,:) - U(end,:)];