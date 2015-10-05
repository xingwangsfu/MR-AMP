classdef CAwgnEstimIn < EstimIn
    % CAwgnEstimIn:  Circular AWGN scalar input estimation function
    
    properties
        % Prior mean and variance
        mean0;  % Mean
        var0;   % Variance
        
        % True indicates to compute output for max-sum
        maxSumVal = false;
    end
    
    methods
        % Constructor
        function obj = CAwgnEstimIn(mean0, var0, maxSumVal)
            obj = obj@EstimIn;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.mean0 = mean0;
                obj.var0 = var0;
                if (nargin >= 3)
                    if (~isempty(maxSumVal))
                        obj.maxSumVal = maxSumVal;
                    end
                end
            end
        end
        
        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = obj.mean0;
            var0  = obj.var0;
            valInit = 0;
        end
        
        % Circular AWGN estimation function
        % Provides the mean and variance of a variable u
        % from an observation v = u + w, w = CN(0,wvar)
        %
        function [umean, uvar, val] = estim(obj, v, wvar)
            % Get prior
            umean0 = obj.mean0;
	    uvar0 = max(eps,obj.var0); % avoid zero variances!
            
            %             % Scale to a vector if they are scalar
            %             if (length(umean0) == 1)
            %                 nv = length(v);
            %                 umean0 = repmat(umean0,nv,1);
            %                 uvar0 = repmat(uvar0,nv,1);
            %             end
            
            % Compute posterior mean and variance
            gain = uvar0./(uvar0+wvar);
            umean = gain.*(v-umean0)+umean0;
            uvar = gain.*wvar;
            
            if (nargout >= 3)
                if ~(obj.maxSumVal) 
                    % Compute the negative KL divergence
                    %   klDivNeg = \sum_i \int p(u|v)*\log( p(u) / p(u|v) )du
		    uvar_over_uvar0 = wvar./(uvar0+wvar);
                    val =  (log(uvar_over_uvar0) + (1-uvar_over_uvar0) ...
                        - abs(umean-umean0).^2./uvar0 );
                else
                    % Evaluate the (log) prior
                    val = -abs(umean-umean0).^2./uvar0;
                end 
            end
            
        end
        
        % Generate random samples
        function x = genRand(obj, outSize)
            if isscalar(outSize)
                x = obj.mean0 +...
                    sqrt(obj.var0/2).*(randn(outSize,1) + 1j*randn(outSize,1));
            else
                x = obj.mean0 +...
                    sqrt(obj.var0/2).*(randn(outSize) + 1j*randn(outSize));
            end
        end
        
        % Computes the likelihood p(y) for y = x + v, v = CN(0,yvar)
        function py = plikey(obj,y,yvar)
            py = exp(-1./((obj.var0+yvar)).*abs(y-obj.mean0).^2);
            py = py./ (pi*(obj.var0+yvar));
        end
        
        % Computes the log-likelihood, log p(y), for y = x + v, where 
        % x = CN(obj.mean0, obj.var0) and v = CN(0, yvar)
        function logpy = loglikey(obj, y, yvar)
            logpy = -( log(pi) + log(obj.var0 + yvar) + ...
                (abs(y - obj.mean0).^2) ./ (obj.var0 + yvar) );
        end
        
    end
    
end

