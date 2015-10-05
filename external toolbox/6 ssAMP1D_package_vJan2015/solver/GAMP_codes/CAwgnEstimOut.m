classdef CAwgnEstimOut < EstimOut
    % CAwgnEstimOut:  CAWGN scalar output estimation function
    %
    % Corresponds to an output channel of the form
    %   y = scale*z + CN(0, wvar)
    
    properties
        % Prior mean and variance
        y;      % Measured output
        wvar;   % Variance
        scale = 1;  % scale factor
        
        % True indicates to compute output for max-sum
        maxSumVal = false;
    end
    
    methods
        % Constructor
        function obj = CAwgnEstimOut(y, wvar, maxSumVal, scale)
            obj = obj@EstimOut;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.y = y;
                obj.wvar = wvar;
                if (nargin >= 3)
                    if (~isempty(maxSumVal))
                        obj.maxSumVal = maxSumVal;
                    end
                end
                if (nargin >= 4)
                    obj.scale = scale;
                end
                
                %Warn user if they set zero variance
                if any(obj.wvar == 0)
                    warning(['Tiny non-zero variances will be used for'...
                        ' computing log likelihoods. May cause problems'...
                        ' with adaptive step size if used.']) %#ok<*WNTAG>
                end
            end
        end
        
        % Size
        function [nz,ncol] = size(obj)
            [nz,ncol] = size(obj.y);
        end
        
        % AWGN estimation function
        % Provides the posterior mean and variance of variable z
        % from an observation y = scale*z + w, z = CN(zmean0,zvar0), w = CN(0,wvar)
        function [zmean, zvar] = estim(obj, zmean0, zvar0)
            
            % Compute posterior mean and variance
            s = obj.scale;
            gain = conj(s)*zvar0./((abs(s)^2)*zvar0 + obj.wvar);
            zmean = gain.*(obj.y-s*zmean0) + zmean0;
            zvar = obj.wvar.*zvar0./((abs(s)^2)*zvar0 + obj.wvar);
            
        end
        
        % Compute log likelihood
        % For sum-product GAMP, compute
        %   E( log p_{Y|Z}(y|z) ) with z = CN(zhat, zvar)
        % For max-sum GAMP compute
        %   log p_{Y|Z}(y|z) @ z = zhat
        function ll = logLike(obj,zhat,zvar)
            
            % Ensure variance is small positive number
            wvar1 = max(1e-20, obj.wvar);
            
            % Get scale
            s = obj.scale;            

            % Compute log-likelihood
            if ~(obj.maxSumVal)
                predErr = (abs(obj.y-s*zhat).^2 + (abs(s)^2)*zvar)./wvar1;
            else
                predErr = (abs(obj.y-s*zhat).^2)./wvar1;
            end
            ll = -(predErr); %return the values without summing
        end
        
        % Compute output cost:
        % For sum-product compute
        %   abs(Axhat-phatfix)^2/(pvar) + log int_z p_{Y|Z}(y|z) CN(z;phatfix, pvar) 
        %   with phatfix such that Axhat=estim(phatfix,pvar).
        % For max-sum GAMP, compute
        %   log p_{Y|Z}(y|z) @ z = Axhat
        function ll = logScale(obj,Axhat,pvar,phat)
                   
            % Ensure variance is small positive number
            wvar1 = max(1e-20, obj.wvar);
            
            %Get the scale
            s = obj.scale;   
           
            % Compute output cost
            if ~(obj.maxSumVal)

                % Compute output cost
                closed_form = true;
                if closed_form
                    
                    %Closed form update
                    ll = -log(abs(s)^2*pvar + wvar1) ...
                        - abs(obj.y - s*Axhat).^2./wvar1 - log(pi); 
                else
                    % Find the fixed-point of phat
                    opt.phat0 = Axhat; % works better than phat
                    opt.alg = 1; % approximate newton's method
                    opt.maxIter = 3; 
                    opt.tol = 1e-4; 
                    opt.stepsize = 1; 
                    opt.regularization = obj.wvar^2;  % works well up to SNR=160dB
                    opt.debug = false;
                    phatfix = estimInvert(obj,Axhat,pvar,opt);

                    % Compute log int_z p_{Y|Z}(y|z) CN(z;phatfix, pvar)
                    ls = -log(pi*(obj.wvar + abs(s)^2*pvar)) ...
                        - abs(obj.y - s*phatfix).^2 ./ (obj.wvar + abs(s)^2*pvar);

                    % Combine to form output cost
                    ll = ls + abs(Axhat - phatfix).^2./pvar;
                end;

            else
                % Output cost is simply the log likelihood
                ll = -abs(obj.y-s*Axhat).^2./wvar1;
            end 
            
        end
        
        function S = numColumns(obj)
            %Return number of columns of Y
            S = size(obj.y,2);
        end
        
        % Generate random samples from p(y|z)
        function y = genRand(obj, z)
            y = sqrt(obj.wvar/2).*randn(size(z)) + ...
                1j*sqrt(obj.wvar/2).*randn(size(z)) + obj.scale.*z;
        end
    end
    
end

