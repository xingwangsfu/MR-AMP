function u = tvdenoise(f,lambda,varargin)
%TVDENOISE  Total variation image denoising.
%   u = TVDENOISE(f,lambda,noise) denoises grayscale, color, or arbitrary
%   multichannel image f using total variation regularization.  Parameter
%   lambda controls the strength of the noise reduction: smaller lambda
%   implies stronger denoising.
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
%   TVDENOISE(...,tol,maxiter) specify the stopping tolerance and the
%   maximum number of iterations.
%
%   TVDENOISE(...,tol,maxiter,plotfun) specifies a plotting callback to
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
%   TVDENOISE(...,tol,maxiter,plotfun,u0) specifies the initial guess u0.  
%   By default, the initialization u0 = f is used.
%
%   See also tvdeconv, tvinpaint, and tvrestore.

% Pascal Getreuer 2009-2010

if nargin < 2
    error('Not enough input arguments.');
end

u = tvrestore(f,lambda,[],[],varargin{:});
