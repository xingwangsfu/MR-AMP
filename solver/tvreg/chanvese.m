function phi = chanvese(f, phi, varargin)
%CHANVESE   Chan-Vese two-phase image segmentation
%   phi = CHANVESE(f,phi0,opt) performs two-phase segmentation on image f 
%   using the Chan-Vese model, where f is the image to be segmented, phi0
%   is the initial level set function, and opt is a struct containing all 
%   or any subset of the following parameters:
%
%      opt.tol       convergence tolerance (default 1e-4)
%      opt.maxiter   maximum number of iterations (default 500)
%      opt.mu        length penalty parameter (default 0.25)
%      opt.nu        area penalty (default 0)
%      opt.lambda1   inside fit penalty (default 1)
%      opt.lambda2   outside fit penalty (default 1)
%      opt.dt        timestep parameter (default 0.5)
%      opt.plotfun   plotting function
%      opt.verbose   show verbose information
%
%   The initial segmentation phi0 should either be phi0=[] for the default
%   or an array such that sign(phi0) indicates the two segments, where 
%   (phi0 >= 0) is the inside and (phi0 < 0) is the outside.
%
%   The plotfun can be used to customize the display of the solution
%   progress.  It should be the name of a function or a function handle.  
%   An example plotfun is
%
%       function myplot(state, iter, delta, phi)
%       switch state
%           case 0  % Running
%               fprintf('  Iter=%4d  Delta=%7.3f\n', iter, delta);
%           case 1  % Converged successfully
%               fprintf('Converged successfully.\n');
%           case 2  % Maximum iterations exceeded
%               fprintf('Maximum number of iterations exceeded.\n');
%       end
%       imagesc(phi >= 0);
%       axis image
%       title(sprintf('Iter=%d  Delta=%.4f', iter, delta));
%       return;
%
%   phi = CHANVESE(f,phi0,tol,maxiter,mu,nu,lambda1,lambda2,dt,plotfun) is
%   an alternative syntax.
%
%   If f is a multichannel image, the segmentation is done using the Chan,
%   Sandberg, Vese vector extension of the Chan-Vese model.

% Pascal Getreuer 2009-2010

% Implementation notes:
%
% This is the pure M-code version of the solver.  It is "functionally
% equivalent" to the MEX function version in chanvese.c, meaning for the 
% same input parameters, they should converge to similar solutions.  
% However, the two versions do not use exactly the same algorithm.
%
% For the semi-implicit update of phi, the MEX version updates the array by
% immediately overwriting coordinate values as they are computed (nonlinear
% Gauss-Seidel).  Here, the M-code version computes all new coordinates
% from the previous coordinate values, and updates them all at once
% (Jacobi).  Jacobi is more friendly for M-code since the operation is
% easily vectorized.


% 
% Input parsing
%

% Check first two arguments f and lambda
if nargin < 1
    error('Not enough input argments.');
elseif (~isnumeric(f) && ~islogical(f)) || ndims(f) > 3
    error('First argument must be a numeric 2D or 3D array.');
end

if nargin < 2 || isempty(phi)
    phi = initdefaultphi(size(f,1), size(f,2));
elseif (~isnumeric(phi) && ~islogical(phi)) || ndims(phi) ~= 2
    error('Second argument must be a numeric 2D array.');
elseif numel(phi) ~= 1 && ~isequal(size(phi),[size(f,1),size(f,2)])
    error('phi0 must have the same number of rows and columns as the input image.');
end

% Set defaults
tol = 1e-4;
maxiter = 500;
mu = 0.25;
nu = 0;
lambda1 = 1;
lambda2 = 1;
dt = 0.5;
plotfun = [];
verbose = 0;

if nargin < 3 
    opt = struct([]);
elseif isstruct(varargin{1})
    if length(varargin) > 1
        error('Too many input arguments');
    end
    
    opt = varargin{1};
else
    if length(varargin) > 8
        error('Too many input arguments');
    end

    opt = {'tol', 'maxiter', 'mu', 'nu', ...
        'lambda1', 'lambda2', 'dt', 'plotfun'};
    opt = cell2struct(varargin, {opt{1:length(varargin)}}, 2);
end

% Unpack opt struct
unpackopt(opt);

% Check scalar parameters
if ~isnumericscalar(mu)
    error('mu should be a numeric scalar.');
elseif ~isnumericscalar(nu)
    error('nu should be a numeric scalar.');
elseif ~isnumericscalar(tol)
    error('tol should be a numeric scalar.');
elseif ~isnumericscalar(maxiter)
    error('maxiter should be a numeric scalar.');
elseif ~isnumericscalar(lambda1)
    error('lambda1 should be a numeric scalar.');
elseif ~isnumericscalar(lambda2)
    error('lambda2 should be a numeric scalar.');
elseif ~isnumericscalar(dt)
    error('dt should be a numeric scalar.');
end

% Check plotfun
if ~isempty(plotfun) && ~ischar(plotfun) ...
        && ~isa(plotfun,'function_handle')
    error('plotfun must be the name of a function or a function handle.');
end

[N1,N2,N3] = size(f);

if verbose > 0
    % Verbose information
    fprintf('f         : [%d x %d x %d]\n', N1, N2, N3);
    fprintf('tol       : %g\n', tol);
    fprintf('max iter  : %d\n', maxiter);
    fprintf('mu        : %g\n', mu);
    fprintf('nu        : %g\n', nu);
    fprintf('lambda1   : %g\n', lambda1);
    fprintf('lambda2   : %g\n', lambda2);
    fprintf('dt        : %g\n', dt);
    fprintf('implement : M-code\n');
end


% 
% Initializations
%

divideeps = 1e-6;

id = [2:N1,N1];
iu = [1,1:N1-1];
ir = [2:N2,N2];
il = [1,1:N2-1];

[c1,c2] = regionaverages(f, phi);

iter = 0;
delta = (tol + (tol == 0))*1000;

if ~isempty(plotfun)
    feval(plotfun, 0, iter, delta, phi);
    drawnow;
end

% 
% Main computation
%

while iter < maxiter
    iter = iter + 1;
    lastphisign = (phi >= 0);
    
    boundary = dt./(1 + phi.*phi);    
    
    PhiX = phi(:,ir) - phi;
    PhiY = ( phi(id,:) - phi(iu,:))/2;
    IDivR = 1./sqrt(divideeps + PhiX.*PhiX + PhiY.*PhiY);
    PhiX = phi - phi(:,il);
    IDivL = 1./sqrt(divideeps + PhiX.*PhiX + PhiY.*PhiY);
    PhiX = (phi(:,ir) - phi(:,il))/2;
    PhiY =  phi(id,:) - phi;
    IDivD = 1./sqrt(divideeps + PhiX.*PhiX + PhiY.*PhiY);
    PhiY = phi - phi(iu,:);
    IDivU = 1./sqrt(divideeps + PhiX.*PhiX + PhiY.*PhiY);

    if N3 == 1
        Dist1 = f - c1;
        Dist1 = Dist1.*Dist1;
        Dist2 = f - c2;
        Dist2 = Dist2.*Dist2;
    else
        Dist1 = 0;
        Dist2 = 0;
        
        for k = 1:N3
            fk = f(:,:,k);
            Temp = fk - c1(k);
            Dist1 = Dist1 + Temp.*Temp;
            Temp = fk - c2(k);
            Dist2 = Dist2 + Temp.*Temp;
        end
    end
    
    % Semi-implicit update
    phi = (phi + boundary.*( ...
            mu*(phi(:,ir).*IDivR + phi(:,il).*IDivL ...
                + phi(id,:).*IDivD + phi(iu,:).*IDivU) ...
            - nu - lambda1*Dist1 + lambda2*Dist2) ) ...
                ./ (1 + boundary.*(mu*(IDivR + IDivL + IDivD + IDivU)));
    
    % Update the region averages c1 and c2
    [c1,c2] = regionaverages(f, phi);
    
    delta = nnz((phi >= 0) ~= lastphisign) / (N1*N2);        
    
    % Convergence check
    if iter >= 2 && delta < tol
        if ~isempty(plotfun)
            feval(plotfun, 1, iter, delta, phi);            
        end
        
        break;
    end
    
    if ~isempty(plotfun)
        feval(plotfun, 0, iter, delta, phi);
        drawnow;
    end
end

if ~isempty(plotfun)
    feval(plotfun, 2, iter, delta, phi);
    drawnow;
end

return;


function phi = initdefaultphi(N1, N2)
% Default initialization for phi
[x,y] = meshgrid(0:N2-1, 0:N1-1);
phi = sin(x*(pi/5)).*sin(y*(pi/5));
return;


function [c1,c2] = regionaverages(f, phi)

[N1,N2,N3] = size(f);
inside = (phi >= 0);

c1 = zeros(1,N3);
c2 = zeros(1,N3);

for k = 1:N3
    fk = f(:,:,k);
    c1(k) = sum(fk(inside));
    c2(k) = sum(fk(~inside));
end

c1 = c1 / max(1,nnz(inside));
c2 = c2 / max(1,nnz(~inside));
return;


function t = isnumericscalar(x)
% Test if variable is a numeric scalar
t = ((isnumeric(x) || islogical(x)) && numel(x) == 1);
return;


function unpackopt(opt)
% Unpack options struct
names = fieldnames(opt);

for k = 1:length(names)
    field = getfield(opt,names{k}); %#ok<GFLD>
    
    if isempty(field)
        continue;
    elseif isempty(strmatch(names{k},{'tol', 'maxiter', 'mu', 'nu', ...
            'lambda1', 'lambda2', 'dt', 'plotfun', 'verbose'}))
        warning(sprintf('Ignoring unknown option ''%s''.', ...
            names{k})); %#ok<WNTAG,SPWRN>
    else
        assignin('caller', names{k}, field);
    end
end
return;
