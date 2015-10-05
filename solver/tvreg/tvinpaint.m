function u = tvinpaint(f,lambda,D,varargin)
%TVINPAINT  Total variation image inpainting.
%   u = TVINPAINT(f,lambda,D,noise) simultaneously denoises and inpaints
%   grayscale, color, or arbitrary multichannel image f using total 
%   variation regularization.  Parameter lambda controls the strength of
%   the noise reduction outside the inpainting region: smaller lambda
%   implies stronger denoising.  D is a logical array specifying the
%   inpainting domain (true = unknown).
%
%   The noise parameter specifies the noise model (case insensitive):
%     'Gaussian' or 'L2'  - (default) The degradation model for additive
%                           white Gaussian noise (AWGN),
%                             f = (exact) + (Gaussian noise).
%     'Laplacian' or 'L1' - The degradation model assumes impulsive noise, 
%                           for example, salt & pepper noise.
%     'Poisson'           - Each pixel is an independent Poisson random
%                           variable with mean equal to the exact value.
%
%   TVINPAINT(...,tol,maxiter) specify the stopping tolerance and the
%   maximum number of iterations.
%
%   TVINPAINT(...,tol,maxiter,plotfun) specifies a plotting callback to
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
%   TVINPAINT(...,tol,maxiter,plotfun,u0) specifies the initial guess u0.  
%   By default, the initialization u0 = f is used.
%
%   See also tvdenoise, tvdeconv, and tvrestore.

% Pascal Getreuer 2009-2010

if nargin < 3
    error('Not enough input arguments.');
end

u = tvrestore(f,lambda,D,[],varargin{:});
