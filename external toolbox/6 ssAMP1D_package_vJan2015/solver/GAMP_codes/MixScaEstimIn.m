classdef MixScaEstimIn < handle
    % MixScaEstimIn:  Scalar estimator class constructed by mixing two other scalar estimators
    %
    % Let X1 and X0 be two random variables.
    % Then MixScaEstimIn is the estimator for a new random variable:
    %   X = { X1 with prob p1
    %       { X0 with prob 1-p1
    properties
        p1;      % prob that estim = estim1
        estim1;  % Base estimator when U=1
        estim0;  % Base estimator when U=1
    end
    
    properties (Hidden)
        LogLikeFlag0 = false;   % estim0 implements loglikey method if true
        LogLikeFlag1 = false;   % estim1 implements loglikey method if true
    end
    
    
    methods
        % Constructor
        function obj = MixScaEstimIn(estim1, p1, estim0)
                obj.p1 = p1;
                obj.estim1 = estim1;
                obj.estim0 = estim0;
        end
        
        % Set method for estim0
        function set.estim0(obj, Estim0)
            % Check to ensure input Estim1 is a valid EstimIn class
            if isa(Estim0, 'EstimIn')
                if ismethod(Estim0, 'loglikey')
                    % Base estimator implements loglikey method
                    obj.LogLikeFlag0 = true;    %#ok<MCSUP>
                end
                obj.estim0 = Estim0;
            else
                error('estim0 must be a valid EstimIn object')
            end
        end
        
        % Set method for estim1
        function set.estim1(obj, Estim1)
            % Check to ensure input Estim1 is a valid EstimIn class
            if isa(Estim1, 'EstimIn')
                if ismethod(Estim1, 'loglikey')
                    % Base estimator implements loglikey method
                    obj.LogLikeFlag1 = true;    %#ok<MCSUP>
                end
                obj.estim1 = Estim1;
            else
                error('estim1 must be a valid EstimIn object')
            end
        end
        
        % Compute prior mean and variance
        function [xhat, xvar, valInit] = estimInit(obj)
            [xhat1, xvar1, valInit1] = obj.estim1.estimInit;
            [xhat0, xvar0, valInit0] = obj.estim0.estimInit;
            xhat = obj.p1.*xhat1 + (1-obj.p1).*xhat0;
            xvar = obj.p1.*(xvar1 + abs(xhat1-xhat).^2) + ...
	    		(1-obj.p1).*(xvar0 + abs(xhat0-xhat).^2);
            valInit = obj.p1.*valInit1 + (1-obj.p1).*valInit0;	% check this!
        end
        
        % Compute posterior outputs
        function [xhat, xvar, klDivNeg, py1] = estim(obj, y, yvar)
            
            % Calculate posterior activity probabilities
            if ~obj.LogLikeFlag0
                % Convert from prob to log-prob domain
                loglike0 = log( obj.estim0.plikey(y, yvar) );
            else
                loglike0 = obj.estim0.loglikey(y, yvar);
            end
            
            if ~obj.LogLikeFlag1
                % Convert from prob to log-prob domain
                loglike1 = log( obj.estim1.plikey(y, yvar) );
            else
                loglike1 = obj.estim1.loglikey(y, yvar);
            end
            
            % Convert log-domain quantities into posterior activity
            % probabilities (i.e., py1 = Pr{X=X1 | y}, py0 = Pr{X=X0 | y})
            exparg = loglike0 - loglike1 + log(1 - obj.p1) - log(obj.p1);
            py1 = (1 + exp(exparg)).^(-1);
            py0 = 1 - py1;
            
            % Compute posterior mean and variance
            [xhat1, xvar1, klDivNeg1] = obj.estim1.estim(y,yvar);
            [xhat0, xvar0, klDivNeg0] = obj.estim0.estim(y,yvar);
            xhat = py1.*xhat1 + py0.*xhat0;
            xvar = py1.*(abs(xhat1-xhat).^2 + xvar1) + ...
                py0.*(abs(xhat0-xhat).^2 + xvar0);
            
            % Compute negative K-L divergence
            if (nargout >= 3)
                klDivNeg = py1.*(klDivNeg1 + log(obj.p1./max(py1,1e-8))) ...
                    + py0.*(klDivNeg0 + log((1-obj.p1)./max(py0,1e-8)));
            end
            
        end
        
        % Generate random samples
        function x = genRand(obj, nx)
            x1 = obj.estim1.genRand(nx);
            x0 = obj.estim0.genRand(nx);
            p = rand(size(x1)) < obj.p1;
            x = x1.*p + x0.*(1-p);
        end
        
        % Get the points in the distribution
        function x0 = getPoints(obj)
            %x0 = [0; obj.estim1.getPoints()];				
            x0 = [obj.estim0.getPoints(); obj.estim1.getPoints()];	% not sure about this!
        end
        
        % Set sparsity level
        function setSparseProb(obj, p1)
            obj.p1 = p1;					% not sure if this needs to be modified!
        end
        
    end
    
end

