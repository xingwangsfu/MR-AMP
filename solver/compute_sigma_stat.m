 function [sigma_stat_final] = compute_sigma_stat (sigma, sigma_e, epision, alpha)

% case 1:  x = 1    0.064
% 1 : z>alpha-1/sigma_stat         % integral region
%  sigma_stat^2*(z-alpha)^2    % the function
% 0 : (-1)*alpha-1/sigma_stat < z < alpha-1/sigma_stat
% 1
% -1 : z < (-1)*alpha-1/sigma_stat
% sigma_stat^2*(z+alpha)^2

% clear all;
sigma_stat_cand = linspace(0.01,1,100);
i = 1;
flag = 0;
while flag ~= 1
  %  sigma_stat = sigma_stat_cand(i);
   % sigma = sqrt(0.2);
   % sigma_e = sqrt(0.2);
  %  alpha = 5;
    %  sigma_stat = 0.63244;
  %  epision = 0.65;
  %   syms sigma_stat;
  sigma_stat = sigma_stat_cand(i);
    syms z ;
    pi = 3.1416;
    func1_1 =   1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
    a1 = int(func1_1, z, alpha-1/sigma_stat,inf);
    func1_2 =  1/sqrt(2*pi)*exp(-z^2/2);
    a2 = int(func1_2, z, (-1)*alpha-1/sigma_stat, alpha-1/sigma_stat);
    func1_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
    a3 = int(func1_3, z, -inf, (-1)*alpha-1/sigma_stat);
    
    % case 0:  x = 0    0.872
    % 1 : z>alpha
    % sigma_stat^2*(z-alpha)^2
    % 0: 0
    
    % -1 : z< (-1)*alpha
    % sigma_stat^2*(z+alpha)^2
    func2_1 =   1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
    b1 = int(func2_1, z, alpha, inf);
    func2_2 = 0;
    b2 = 0;
    func2_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
    b3 = int(func2_3, z, -inf, (-1)*alpha);
    
    % case -1: x = -1    0.064
    % 1: z>alpha+1/sigma_stat
    % sigma_stat^2*(z-alpha)^2
    
    % 0: (-1)*alpha+1/sigma_stat < z < alpha+1/sigma_stat
    % 1
    % -1: z < (-1)*alpha+1/sigma_stat
    % % sigma_stat^2*(z+alpha)^2
    
    func3_1 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
    c1 = int(func3_1, z, alpha+1/sigma_stat, inf);
    func3_2 =  1/sqrt(2*pi)*exp(-z^2/2);
    c2 = int(func3_2, z, (-1)*alpha+1/sigma_stat, alpha+1/sigma_stat);
    func3_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
    c3 = int(func3_3, z, -inf, (-1)*alpha+1/sigma_stat);
    
    fMSE = (a1+a2+a3)*0.064+(b1+b2+b3)*0.872+(c1+c2+c3)*0.064;
    eq1(i) = sigma_stat^2-sigma_e^2*(sigma^2+fMSE/epision)/(sigma_e^2+sigma^2+fMSE/epision);
 %   sigma_stat_final  = solve(eq1,sigma_stat);
    if i>=2 
        if double(eq1(i-1))<=0 && double(eq1(i))>=0 
        flag = 1;
        end
    end
    equ_rhs(i) = sigma_e^2*(sigma^2+fMSE/epision)/(sigma_e^2+sigma^2+fMSE/epision);
    i = i+1;
end
[v, I] = min(abs(double(eq1)));
sigma_stat_final = sigma_stat_cand(I);
% equ_rhs = sigma_e^2*(sigma^2+fMSE/epision)/(sigma_e^2+sigma^2+fMSE/epision);
%  ppp = solve(eq1,sigma_stat);