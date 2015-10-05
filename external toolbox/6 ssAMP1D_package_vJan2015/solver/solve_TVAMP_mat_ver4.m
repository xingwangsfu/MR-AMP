%% /*======================================================================
%  * solve_TVAMP_mat_ver4.m: solver of TV-AMP algorithm using EFLA
%  algorithm
%  * 
%  *  This is an AMP algorithm to solver TV minimizatio problem
%  *  We refer to Donoho et al.' paper published in IEEE inform.theory 2013
%  *  For threshoding function in TV-AMP, we have used flsa solver in the SLEP
%  *   package introduced in the paper, 
%  *   [J. Liu, L. Yuan, and J. Ye. ¡°An efficient algorithm for a class
%  *   of fused lasso problems,¡± proc of ACM SIGKDD Conference on
%  *   Knowledge Discovery and Data Mining, 2010.]
%  *   And the code of flsa solver is available at http://www.public.asu.edu/?jye02/Software/SLEP
%  * 
%  * Input parameters (A,y,maxiter,EMopt,DeltaN)
%  * 1) A (M by N): sparse measurement matrix 
%  * 2) y (M by 1): given measurement vector with size M
%  * 3) maxiter   : the number of iterations of the algorithm. 
%  * 4) X0         : true value of signal is used to find minimax
%  *                thresholding paramter lambda
%  * 5) tol    : stopping tolerance - ssAMP will be terminated when the nomarlized MSE <= tol 
%  * 
%  * Output parameters  [meanX,varX]
%  * 1) x_hat (N by 1) : the mean of the marginal posterior density
%  * of the signal, which is an estimate of the signal. 
%  * 2) tau (N by maxiter) : the variance of the marginal posterior
%  * density of the signal, which shows the convergence of the BP
%  * iteration.
%  *
%  * copyright@ Jaewook Kang with Gwangju Institute of Science and
%  * Technology 2014,Jan,  
%  * feedback: jwkkang@gist.ac.kr
%  * final update Mar. 25th, 2014
%%  *=======================================================================*/
function [x,tau]=solve_TVAMP_mat_ver4(A,y,maxiter,X0,tol)



[M,N]=size(A);
x=zeros(N,1);
prev_x=zeros(N,1);

prev_z=zeros(M,1);
z=zeros(M,1);
tau=zeros(maxiter,1);
%initialization
prev_z=y;

% parameter for tvdip solver
maxiterTV=1000;
% lratio = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10];
 lratio = [1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 1e-1 0.5 1 5 10];
%lratio=[10 1 5.^-(1:10)];
% % lratio=[10 1 2.^-(1:20)];
b0=zeros(N-1,1);

for t=2:maxiter+1
    rho=A'*prev_z+prev_x;
     %------------------- thresholding via flsa solver--------------------%
     prev_x=x;
     xMSE=1;
     for tt=1:length(lratio )
        lambda=lratio(tt);
        [temp_x, b, infor]=flsa(rho, b0,0, lambda, N,maxiterTV, 1e-5, 1, 6);
        temp_xMSE=norm(X0-temp_x)^2/norm(X0)^2;
        
        if xMSE > temp_xMSE
            xMSE=temp_xMSE;
            b0=b;
            x=temp_x;
        end
        
     end 
    % termination condition 
    if (norm(prev_x - x)^2/norm(prev_x)^2 <tol) &&(t>20)
        break;
    end
    %--------------------------------------------------------------------%
    % Onsager term calculation 
    Onsager = 1/M* (length(nonzeros(diff(prev_x))));  
%         Onsager = 1/M* (length(nonzeros(diff(prev_x))))-1;  for small M/N
     % residual message update
   %  z = 0.3*prev_z +  0.7*(y - A*x +  prev_z * Onsager) ;
%      z = 0.4*prev_z +  0.6*(y - A*x +  prev_z * Onsager) ;
   z = 0.5*prev_z +  0.5*(y - A*x +  prev_z * Onsager) ;
     prev_x=x;
     prev_z=z;
end

% disp(sprintf('Number of TV-AMP iterations = %d',t-1));

end


    