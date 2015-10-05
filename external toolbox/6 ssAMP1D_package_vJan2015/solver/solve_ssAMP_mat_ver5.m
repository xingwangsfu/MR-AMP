%% --------------------------------------------------------------
%  filename: solve_ssAMP_mat_ver5.m -ss AMP algorithm (m-code
%  version)
% This algorithm is AMP algorithm for the total variation problem.
% The message exchange between m-factor and variable nodes are done in the
% way of AMP.
% The message exchange between s-factor and variable nodes are done in the
% way of parameter MP
%-------------------------------------------------------------------------
% In version 5, complexity reduction is performed from the version 3. 
%-------------------------------------------------------------------------
% <Input parameters (H,y,Delta,sigma,q,maxiter)>
% * 1) H (M by N): measurement matrix
% * 2) y (M by 1): The given measurement vector with size M
% * 3) sigma:  signal variance from the prior knowledge
% * 4) q    : a rate of piecewise-constantness 
% * 5) Delta:  The variance of additive noise
% * 6) maxiter: Allowdable maximum number of AMP iterations
% * 7) tol    : stopping tolerance - ssAMP will be terminated when the nomarlized MSE <= tol 

% <Output parameters (mu,sigmapow,end_iter,theta)>
% * 1) mu: signal estimate from the AMP iteration 
% * 2) sigmapow: the variance of marginal posterior pdf of estimate signals 
% * 3) end_iter: the number of iterations performed in recovery.
% * 4) theta   : the variance of noise-like (noise + estimation error)
%                component at the end of the AMP iteration
% * 5) MSE_iter: MSE of the ssAMP estimate over iteration
% updated at 2014 Feb. Jaewook Kang
%% --------------------------------------------------------------
function [mu,sigmapow,end_iter,theta]= solve_ssAMP_mat_ver5(H,y,sigma0,q,Delta,maxiter,tol)

if nargin < 6
    maxiter=500;
end
if nargin < 5
    Delta = 0.01;
end

sigmapow0=sigma0^2;
% parameter declear for message passing iteration 
[M,N]=size(H);

% the data space declear for messages
Z         =zeros(N,1);
mu         = zeros(N,1);  % signal estimate
prev_mu    =zeros(N,1);
sigmapow  = ones(N,1)*sigmapow0;   % estimate variance


theta = zeros(1, maxiter+1); % 

A=zeros(N,4);
B=zeros(N,4);
C=zeros(N,4);
F=zeros(N,4);

theta_plus_B=zeros(N,4);
normpdf_matrix=zeros(N,4);

theta_plus_sigmapow=zeros(N,2);
normpdfRL_matrix=zeros(N,2);


% messages for left-right passing
Z_left = zeros(N,1);
Z_right = zeros(N,1);
mu_left  = zeros(N,1);
mu_right = zeros(N,1);
sigmapow_left  = sigmapow;
sigmapow_right = sigmapow;

mu_left_prev  = zeros(N,1);
mu_right_prev = zeros(N,1);
sigmapow_left_prev  = sigmapow;
sigmapow_right_prev = sigmapow;

rho_shift=zeros(N,1);
rho      =zeros(N,1);

% initialization
theta(1) = Delta + sum(sigmapow(1))/M;
r        = y; % residual vector

for t=2:maxiter+1
    
   
   %disp(sprintf('//--------------------------The %d -th iteration----------------------//',t-1))
    %% Left/Right message update
    rho      = H'*r+mu;
    theta(t) = Delta + sum(sigmapow)/M ;   
    
    % left-message passing update (i.e. the message passing is toward left)
    mu_left_prev(1:N-1) = mu_left(2:N);    sigmapow_left_prev(1:N-1) = sigmapow_left(2:N);
    mu_right_prev(2:N)  = mu_right(1:N-1); sigmapow_right_prev(2:N)  = sigmapow_right(1:N-1); 
    rho_shift(1:N-1)    = rho(2:N);   rho_shift(N)=0; 

    %-----------------------------------------------------------------------------------------------------------------------------------------------------
    theta_plus_sigmapow = [theta(t)+sigmapow_left_prev, theta(t)+sigmapow_left_prev+sigmapow0];
    normpdfRL_matrix    = [exp(-0.5 * (rho_shift - mu_left_prev).^2./theta_plus_sigmapow(:,1)) ./ sqrt(2*pi*theta_plus_sigmapow(:,1)),...
                           exp(-0.5 * (rho_shift - mu_left_prev).^2./theta_plus_sigmapow(:,2)) ./ sqrt(2*pi*theta_plus_sigmapow(:,2))]*sqrt(theta(t));
    %-----------------------------------------------------------------------------------------------------------------------------------------------------    
    
    Z_left         = zetaRLfunc_AMP_ver2(normpdfRL_matrix,q);
    mu_left       = etaRLfunc_AMP_ver2 (rho_shift,theta(t),mu_left_prev,sigmapow_left_prev,normpdfRL_matrix,theta_plus_sigmapow,q,sigmapow0,Z_left);
    mu_left(N)=0; 
    
    sigmapow_left = gammaRLfunc_AMP_ver2 (rho_shift,theta(t),mu_left_prev,sigmapow_left_prev,normpdfRL_matrix,theta_plus_sigmapow,q,sigmapow0,Z_left,mu_left);
    sigmapow_left(N)=sigmapow0;
    
     % right-message passing update (i.e. the message passing is toward right)
    rho_shift(2:N)      = rho(1:N-1);    rho_shift(1)=0;  
    %-----------------------------------------------------------------------------------------------------------------------------------------------------
    theta_plus_sigmapow = [theta(t)+sigmapow_right_prev, theta(t)+sigmapow_right_prev+sigmapow0];
    normpdfRL_matrix    = [exp(-0.5 * (rho_shift - mu_right_prev).^2./theta_plus_sigmapow(:,1)) ./ sqrt(2*pi*theta_plus_sigmapow(:,1)),...
                           exp(-0.5 * (rho_shift - mu_right_prev).^2./theta_plus_sigmapow(:,2)) ./ sqrt(2*pi*theta_plus_sigmapow(:,2))]*sqrt(theta(t));
    %-----------------------------------------------------------------------------------------------------------------------------------------------------                   
    Z_right        = zetaRLfunc_AMP_ver2(normpdfRL_matrix,q);
    mu_right       = etaRLfunc_AMP_ver2 (rho_shift,theta(t),mu_right_prev,sigmapow_right_prev,normpdfRL_matrix,theta_plus_sigmapow,q,sigmapow0,Z_right);
    mu_right(1)=0;
    
    sigmapow_right = gammaRLfunc_AMP_ver2 (rho_shift,theta(t),mu_right_prev,sigmapow_right_prev,normpdfRL_matrix,theta_plus_sigmapow,q,sigmapow0,Z_right,mu_right);
    sigmapow_right(1)=sigmapow0;
  

        
    %% A,B,C matrices construction
    var_sum1 = sigmapow_left+sigmapow_right;
    var_sum2 = var_sum1 + sigmapow0;
    var_sum4 = var_sum2 + sigmapow0;
    
    
    var_sumL=  sigmapow_left + sigmapow0;
    var_sumR=  sigmapow_right + sigmapow0;
    
    muR_mul_var_sumL = mu_right.* var_sumL;
    muR_mul_varL     = mu_right.* sigmapow_left;
    
    muL_mul_var_sumR = mu_left.*  var_sumR;
    muL_mul_varR     = mu_left.*  sigmapow_right;
    mean_diff_sqr=(mu_left - mu_right).^2;
    
    A =[(muL_mul_varR       + muR_mul_varL) ./var_sum1 ,...
        (muL_mul_var_sumR   + muR_mul_varL) ./var_sum2 ,...
        (muL_mul_varR       + muR_mul_var_sumL)./var_sum2 ,...
        (muL_mul_var_sumR   + muR_mul_var_sumL)./var_sum4 ];
    
    B = [sigmapow_left  .* sigmapow_right  ./var_sum1 ,...
         sigmapow_left  .* var_sumR        ./var_sum2 ,...
         var_sumL       .* sigmapow_right  ./var_sum2 ,...
         var_sumL       .* var_sumR        ./var_sum4 ];
    
    c23= exp(- mean_diff_sqr/2./var_sum2 )./sqrt(var_sum2 );
    C = [ exp(- mean_diff_sqr/2./var_sum1 )./sqrt(var_sum1 ),...
          c23                                               ,...
          c23                                               ,...
          exp(- mean_diff_sqr/2./var_sum4 )./sqrt(var_sum4 )]/sqrt(2*pi);
    clear c23 
%------------------------------------------------------------------------------------------------------------------
    F = [ theta(t).*A(:,1) + rho.*B(:,1),...
          theta(t).*A(:,2) + rho.*B(:,2),...
          theta(t).*A(:,3) + rho.*B(:,3),...
          theta(t).*A(:,4) + rho.*B(:,4)];  
%       
   theta_plus_B=[theta(t)+B(:,1),...
                 theta(t)+B(:,2),...
                 theta(t)+B(:,3),...
                 theta(t)+B(:,4)];

    normpdf_matrix=[exp(-0.5 * (rho - A(:,1)).^2./theta_plus_B(:,1)) ./ sqrt(2*pi*theta_plus_B(:,1)),...
                    exp(-0.5 * (rho - A(:,2)).^2./theta_plus_B(:,2)) ./ sqrt(2*pi*theta_plus_B(:,2)),...
                    exp(-0.5 * (rho - A(:,3)).^2./theta_plus_B(:,3)) ./ sqrt(2*pi*theta_plus_B(:,3)),...
                    exp(-0.5 * (rho - A(:,4)).^2./theta_plus_B(:,4)) ./ sqrt(2*pi*theta_plus_B(:,4))]*sqrt(theta(t));              
%---------------------------------------------------------------------------------------------------------------------
    %% signal to measurements message update
    prev_mu=mu;

    Z         = zetafunc_AMP_ver2(normpdf_matrix,C,q);
    mu        = etafunc_AMP_ver2(normpdf_matrix,theta_plus_B,C,F,q,Z);
    
    % variance calculation 
   sigmapow =  gammafunc_AMP_ver2(theta(t),normpdf_matrix,theta_plus_B,B,C,F,q,Z,mu);
   

    % loop termination condition check 
    if (norm(prev_mu - mu)^2/norm(prev_mu)^2 <tol) 
        break;
    end
    


     %%  measurement to signal message update 
     r = 0.05*r + 0.95*( y - H*mu + r.*sum(DIV_etafunc_AMP_ver2(rho,normpdf_matrix,theta_plus_B,A,B,C,F,q,Z,mu))/M);
%          r = 0.5*r + 0.5*( y - H*mu + r.*sum(DIV_etafunc_AMP_ver2(rho,normpdf_matrix,theta_plus_B,A,B,C,F,q,Z,mu))/M);
end
end_iter=t-1;
disp(sprintf('Number of ssAMP iterations = %d', end_iter));
end


