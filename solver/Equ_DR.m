function [equ_dr] = Equ_DR(sigma_stat, alpha, h_posnega, epsilon)

% to get the equilibrium detection rate
syms z;
pi = 3.1415926;
func =  1/sqrt(2*pi)*exp(-z^2/2);
a1_dr  = int(func, z, (-1)*alpha-h_posnega/sigma_stat, alpha-h_posnega/sigma_stat);
a1_conj = 1-vpa(a1_dr);

% 0
a2_dr = int(func, z, (-1)*alpha, alpha);
a2_conj = 1-vpa(a2_dr);

% -1
a3_dr = int(func, z, (-1)*alpha+h_posnega/sigma_stat, alpha+h_posnega/sigma_stat);
a3_conj = 1-vpa(a3_dr);

% 
equ_dr = epsilon/2*a1_conj+(1-epsilon)*a2_conj+epsilon/2*a3_conj;
