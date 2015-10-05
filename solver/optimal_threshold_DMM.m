function [beta,alpha_num_pos,rho_num] = optimal_threshold_DMM(delta)

% this function is to compute the optiam threhsolding parameter given the
% sparsity epsilon
syms beta;
psi = exp(-1*(beta^2) /2)/sqrt(2*pi);

Phi = int(psi, -inf, -1*beta);

rho = (1-2/delta*((1+beta^2)*Phi-beta*psi))/(1+beta^2-2*((1+beta^2)*Phi-beta*psi));

diff_f = diff(rho);

alpha_num = solve(diff_f, beta);

alpha_num_pos =  abs(double(alpha_num));
rho_num = subs(rho,beta,abs(double(alpha_num)));