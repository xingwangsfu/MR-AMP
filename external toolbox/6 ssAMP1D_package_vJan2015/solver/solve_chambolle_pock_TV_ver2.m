%------------------------------------------------------------------------%
% filename :  solve_chambolle_pock_TV.m
%  This solver aim to solver the 1-dimensional TV regularization problem (unconstraint). 
%
%         min_x ||Ax - b ||_2^2 + lambda || x ||_TV
%
%  We refer to the papaer 
%  Sidhy et al, "Convex optimization problem portotyping for image
%  reconstruction in computed tomography with the Chambolle-Pock
%  algorithm", 2012 physics in medicine and biology.
%  * ----------------------------------------------------------------------
%   In the version 2, the solver does not preserves value of the variables 
%   at the previous iteration. 
%  * ----------------------------------------------------------------------
%  *Input parameters
%  - A (M by N matrix)        : measurement matrix 
%  - b (M by 1 vector)        : measurement vector
%  - x0 (N by 1 vector)       : initial position of x
%  - lambda (scalar )         : The TV regularization parameter
%  - epsilon (scalar)         : termination condition (stopping tolerance )
%  - maxiter (scalar)         : the maxinum number of iterations
%  - L                         : maximum signgular value 
%
%  *output parameters
% - sol_x                     : solution of the optimization 
% -
% copyright@ Jaewook Kang, GIST-CSNL 2014, May 15th
%-------------------------------------------------------------------------%
function [est_x]= solve_chambolle_pock_TV_ver2(A,b,x0,lambda,epsilon,maxiter,L)

[M,N]=size(A);
D=sparse(N,N);

curr_x=zeros(N,1);
curr_p=zeros(M,1);
curr_q=zeros(N,1);
curr_x_bar=zeros(N,1);

next_x=zeros(N,1);
next_p=zeros(M,1);
next_q=zeros(N,1);
next_x_bar=zeros(N,1);
lambda_mul_one_vectorN=ones(N,1)*lambda;
% difference matrix (also called TV norm matrix )
I = speye(N-1);
P = spalloc(N-1,1,0);
D(1:N-1,:)=( [P I] - [I P] );
clear I P



% initialization
curr_x_bar=x0;
tau=1/L;
sigma=1/L;
theta=1;

%--for efficient computation ---%
Dt=D';
tauAt=tau*A';
sigmaD=sigma*D;
%---------------------
% iteration start
for n=1:maxiter
   next_p    = (curr_p+sigma*(A*curr_x_bar-b))/(1+sigma);
   temp1=curr_q+sigmaD*curr_x_bar ;
   next_q    = lambda* temp1 ./ max(lambda_mul_one_vectorN, abs( temp1));
   next_x     = curr_x-tauAt*next_p - tau*Dt*next_q;
   next_x_bar = next_x + theta*(next_x-curr_x);
  
   
   % stopping criterion
   if norm(curr_x-next_x,2)^2/norm(curr_x,2)^2  <epsilon
%        disp(sprintf('Number of Chambolle-Pock iterations = %d\n',n));
       break;
   end
   
   curr_x=next_x;
   curr_x_bar=next_x_bar;
   curr_p=next_p;
   curr_q=next_q; 
end

est_x=next_x;
disp(sprintf('Number of Chambolle-Pock iterations = %d',n));
end
