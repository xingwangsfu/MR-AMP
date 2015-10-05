clear all;

n = 256*256;
delta = 0.2;
m = floor(delta*n);

x = randn(n,1);
y = phi_fp(x,m,n);

y_1 = phi_fp(x,m,n);