classdef SparseScaEstim < handle
    % SparseScaEstim:  Scalar estimator class with sparsity
    %
    % Let baseEstim be the estimator for a base random variable X1.
    % Then SparseScaEstim is the estimator for a new random variable:
    %   X = X1 with prob p1
    %     = x0 with prob 1-p1
    % where x0 is a constant (default = 0)
    properties
        p1;      % Sparsity
        estim1;  % Base estimator when U=1
        x0;      % Null value
    end
    
    properties (Hidden)
        LogLikeFlag = false;    % Indicates if estim1 has a loglikey method
    end
    
    methods
        % Constructor
        function obj = SparseScaEstim(estim1, p1, x0)
            if (nargin < 3)
                x0 = 0;
            end
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.p1 = p1;
                obj.estim1 = estim1;
                obj.x0 = x0;
            end
        end
        
        % Set method for estim1
        function set.estim1(obj, Estim1)
            % Check to ensure input Estim1 is a valid EstimIn class
            if isa(Estim1, 'EstimIn')
                if ismethod(Estim1, 'loglikey')
                    % Base estimator implements loglikey method
                    obj.LogLikeFlag = true; %#ok<MCSUP>
                end
                obj.estim1 = Estim1;
            else
                error('estim1 must be a valid EstimIn object')
            end
        end
        
        % Compute prior mean and variance
        function [xhat, xvar, valInit] = estimInit(obj)
            [xhat1, xvar1, valInit1] = obj.estim1.estimInit;
            xhat = obj.p1.*xhat1 + (1-obj.p1).*obj.x0;
            xvar = obj.p1.*(abs(xhat1).^2 + xvar1) + (1-obj.p1).*(abs(obj.x0).^2) - ...
                abs(xhat).^2;
            valInit = obj.p1.*valInit1;
        end
        
        % Compute posterior mean and variance from Gaussian estimate
        function [xhat, xvar, klDivNeg, py1] = estim(obj, y, yvar)
            
            % Compute the activity probabilities
            if ~obj.LogLikeFlag
                % This EstimIn class does not implement the loglikey
                % method, thus convert from probability domain to log-prob
                % domain
                
                % Get log-likelihood of y for U=1 and U=0
                loglike1 = log( obj.estim1.plikey(y, yvar) );

            else
                % This EstimIn class implements the loglikey method, thus
                % work in the log-domain
                
                % Get log-likelihood of y for U=1 and U=0
                loglike1 = obj.estim1.loglikey(y, yvar);
            end
            
            %Handle real and complex cases separately
            yvar(yvar < eps) = eps;     % for numerical stability...
            if isreal(y)
                loglike0 = (-1/2)*( log(2*pi) + log(yvar) + ...
                    ((y - obj.x0).^2)./yvar );
            else
                loglike0 = -( log(pi) + log(yvar) + ...
                    (abs(y - obj.x0).^2)./yvar );
            end
            
            % Convert log-domain quantities into posterior activity
            % probabilities (i.e., py1 = Pr{X=X1 | y}, py0 = Pr{X=x0 | y})
            exparg = loglike0 - loglike1 + log(1 - obj.p1) - log(obj.p1);
            py1 = (1 + exp(exparg)).^(-1);
            py0 = 1 - py1;
            
            % Compute mean and variance
            if (nargout >= 3)
                [xhat1, xvar1, klDivNeg1] = obj.estim1.estim(y, yvar);
            else
                [xhat1, xvar1] = obj.estim1.estim(y, yvar);
            end
            xhat = py1.*xhat1 + py0.*obj.x0;
            xvar = py1.*(abs(xhat1).^2 + xvar1) + py0.*(abs(obj.x0).^2)...
                - abs(xhat).^2;
            
            % Compute negative K-L divergence
            if (nargout >= 3)
                klDivNeg = py1.*klDivNeg1 + py1.*log(max(1e-8,obj.p1)./max(py1,1e-8)) ...
                    + py0.*log(max(1e-8,(1-obj.p1))./max(py0,1e-8));
            end
            
        end
        
        % Generate random samples
        function x = genRand(obj, nx)
            x1 = obj.estim1.genRand(nx);
            p = rand(size(x1)) < obj.p1;
            x = x1.*p + (1-p).*obj.x0;
        end
        
        % Get the points in the distribution
        function x0 = getPoints(obj)
            x0 = [0; obj.estim1.getPoints()];
        end
        
        % Set sparsity level
        function setSparseProb(obj, p1)
            obj.p1 = p1;
        end
        
    end
    
end

