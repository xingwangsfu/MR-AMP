function u = tvdeconv(f,lambda,K,varargin)
%TVDECONV  Total variation image deconvolution.
%   u = TVDECONV(f,lambda,K,noise) deblurs the input image f such that u
%   is approximately K*f.  The parameter lambda balances between
%   deblurring accuracy and noise reduction, where smaller lambda implies
%   stronger noise reduction (but at the cost of deblurring accuracy).  The
%   method requires that K has nonzero sum, sum(K(:)) ~= 0.
%
%   The noise parameter specifies the noise model (case insensitive):
%     'Gaussian' or 'L2'  - (default) the degradation model for additive
%                           white Gaussian noise (AWGN),
%                             f = (exact) + (Gaussian noise)
%     'Laplace' or 'L1'   - additive Laplacian noise, this model is
%                           effective for salt&pepper noise
%     'Poisson'           - each pixel is an independent Poisson random
%                           variable with mean equal to the exact value
%
%   TVDECONV(...,tol,maxiter) specify the stopping tolerance and the
%   maximum number of iterations.
%
%   TVDECONV(...,tol,maxiter,plotfun) specifies a plotting callback to
%   customize the display of the solution progress.  plotfun should be the
%   name of a function or a function handle.  An example plotfun is
%
%       function myplot(state, iter, delta, u)
%       switch state
%           case 0  % Running
%               fprintf('  Iter=%4d  Delta=%7.3f\n', iter, delta);
%           case 1  % Converged successfully
%               fprintf('Converged successfully.\n');
%           case 2  % Maximum iterations exceeded
%               fprintf('Maximum number of iterations exceeded.\n');
%       end
%       image(u*255);
%       axis image
%       title(sprintf('Iter=%d  Delta=%.4f', iter, delta));
%       return;
%
%   TVDECONV(...,tol,maxiter,plotfun,u0) specifies the initial guess u0.  
%   By default, the initialization u0 = f is used.
%
%   See also tvdenoise, tvinpaint, and tvrestore.

% Pascal Getreuer 2009-2010

if nargin < 3
    error('Not enough input arguments.');
end

u = tvrestore(f,lambda,[],K,varargin{:});
