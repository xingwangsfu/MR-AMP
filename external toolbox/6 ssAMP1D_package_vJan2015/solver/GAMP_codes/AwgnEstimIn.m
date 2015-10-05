classdef AwgnEstimIn < EstimIn
    % AwgnEstimIn:  AWGN scalar input estimation function
    
    properties 
        % Prior mean and variance
        mean0;  % Mean 
        var0;   % Variance 

        % True indicates to compute output for max-sum
        maxSumVal = false;
    end
    
    methods
        % Constructor
        function obj = AwgnEstimIn(mean0, var0, maxSumVal)
            obj = obj@EstimIn;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.mean0 = mean0;
                obj.var0 = var0;
                if (nargin >= 3)
                    if (~isempty(maxSumVal))
                        obj.maxSumVal = maxSumVal;
                    end
                end
                
                % warn user about inputs
                if any(~isreal(mean0(:))),
                    error('First argument of AwgnEstimIn must be real-valued');
                end;
                if any((var0(:)<0))||any(~isreal(var0(:))),
                    error('Second argument of AwgnEstimIn must be non-negative');
                end;
            end
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = obj.mean0;
            var0  = obj.var0;
            valInit = 0;
        end

        % AWGN estimation function
        % Provides the mean and variance of a variable u 
        % from a observation real(v) = u + w, w = N(0,wvar)
        %
        function [umean, uvar, val] = estim(obj, v, wvar)
            % Get prior
            umean0 = obj.mean0;
	    uvar0 = max(eps,obj.var0); % avoid zero variances!
            
            % Compute posterior mean and variance
            gain = uvar0./(uvar0+wvar);
            umean = gain.*(real(v)-umean0)+umean0;
            uvar = gain.*wvar;
            
            if (nargout >= 3)                            
                if ~(obj.maxSumVal)
                    % Compute the negative KL divergence            
                    %   klDivNeg = \sum_i \int p(u|v)*\log( p(u) / p(u|v) )du
		    uvar_over_uvar0 = wvar./(uvar0+wvar);
                    val = 0.5* (log(uvar_over_uvar0) + (1-uvar_over_uvar0) ...
                        - (umean-umean0).^2./uvar0 );
                else
                    % Evaluate the (log) prior
                    val = -0.5* (umean-umean0).^2./uvar0;
                end
            end

        end
        
        % Generate random samples
        function x = genRand(obj, outSize)
            if isscalar(outSize)
                x = sqrt(obj.var0).*randn(outSize,1) + obj.mean0;
            else
                x = sqrt(obj.var0).*randn(outSize) + obj.mean0;
            end
        end
        
        % Computes the likelihood p(y) for real(y) = x + v, v = N(0,yvar)
        function py = plikey(obj,y,yvar)
            py = exp(-1./(2*(obj.var0+yvar)).*(real(y)-obj.mean0).^2);
            py = py./ sqrt(2*pi*(obj.var0+yvar));
        end
        
        % Computes the log-likelihood, log p(y), for real(y) = x + v, where 
        % x = N(obj.mean0, obj.var0) and v = N(0, yvar)
        function logpy = loglikey(obj, y, yvar)
            logpy = (-0.5)*( log(2*pi) + log(obj.var0 + yvar) + ...
                ((real(y) - obj.mean0).^2) ./ (obj.var0 + yvar) );
        end

    end
    
end

