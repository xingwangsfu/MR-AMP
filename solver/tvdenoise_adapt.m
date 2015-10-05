function u = tvdenoise_adapt(f,sigma)

N = length(f);
n = sqrt(N);
f = reshape(f,n,n);
f = f/255;
sigma = sigma/255;
% initial lambda

lambda = 0.7079/sigma + 0.6849/(sigma^2);

for i = 1:5
    u = tvdenoise(f,lambda);
    tmp = f(:)-u(:);
    tmp_lambda(i) = lambda;
    lambda = lambda*sqrt(sum(tmp.^2)/length(u(:)))/sigma;
    if abs(tmp_lambda(i)-lambda)/tmp_lambda(i) < 1e-2
        break;
    end
    % norm(i) = sum((u(:)-X(:)).^2)/(length(u(:)));
end
u = u(:)*255;
