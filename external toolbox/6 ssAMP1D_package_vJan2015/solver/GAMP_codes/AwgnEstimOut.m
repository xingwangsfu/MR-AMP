classdef AwgnEstimOut < EstimOut
    % AwgnEstimOut:  AWGN scalar output estimation function
    %
    % Corresponds to an output channel of the form
    %   y = scale*z + N(0, wvar)
    
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
        function obj = AwgnEstimOut(y, wvar, maxSumVal, scale)
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
                    if (scale <= 0),
                        error('Fourth argument of AwgnEstimOut must be positive');
                    end;
                end
                
                %Warn user about inputs
                if any(~isreal(obj.y)),
                    error('First argument of AwgnEstimOut must be real-valued.  Did you mean to use CAwgnEstimOut instead?');
                    % if we really want to handle real-valued noise and complex-valued y
                    % (and thus complex z), then we need to modify this file!
                end;
                if any(any((obj.wvar<0)))||any(any(~isreal(obj.wvar))),
                    error('Second argument of AwgnEstimOut must be non-negative');
                end;
                if any(obj.wvar==0)
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
        % Provides the posterior mean and variance of _real_ variable z
        % from an observation real(y) = scale*z + w
        % where z = N(zmean0,zvar0) and w = N(0,wvar)
        function [zmean, zvar] = estim(obj, zmean0, zvar0)
            
            % Compute posterior mean and variance
            s = obj.scale;
            gain = conj(s)*zvar0./((s^2)*zvar0 + obj.wvar);
            zmean = gain.*(obj.y-s*real(zmean0)) + real(zmean0);
            zvar = obj.wvar.*zvar0./((s^2)*zvar0 + obj.wvar);
            
        end
        
        % Compute output cost:
        % For sum-product compute
        %   E_Z( log p_{Y|Z}(y|z) ) with Z ~ N(zhat, zvar)
        % For max-sum GAMP, compute
        %   log p_{Y|Z}(y|z) @ z = zhat
        function ll = logLike(obj,zhat,zvar)
            
            % Ensure variance is small positive number
            wvar1 = max(eps, obj.wvar);
            
            % Get scale
            s = obj.scale;
            
            % Compute log-likelihood
            if ~(obj.maxSumVal)
                predErr = ((obj.y-s*real(zhat)).^2 + (s^2)*zvar)./wvar1;
            else
                predErr = ((obj.y-s*real(zhat)).^2)./wvar1;
            end
            ll = -0.5*(predErr); %return the values without summing
        end
        
        % Compute output cost:
        % For sum-product compute
        %   (Axhat-phatfix)^2/(2*pvar) + log int_z p_{Y|Z}(y|z) N(z;phatfix, pvar) 
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
                    %ptil = (pvar./wvar1+1).*Axhat - (pvar./wvar1).*obj.y;
                    %Old closed form update without scale
%                     ll = -0.5*( log(pvar+wvar1) + (obj.y-real(Axhat)).^2./wvar1 ...
%                         + log(2*pi) );
                    
                    % Compute in closed form
                    ll = -0.5*( log(s^2*pvar + wvar1) ...
                        + (obj.y - s*real(Axhat)).^2./wvar1 + log(2*pi)); 
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

                    % Compute log int_z p_{Y|Z}(y|z) N(z;phatfix, pvar)
                    ls = -0.5*log(2*pi*(obj.wvar + s^2*pvar)) ...
                        - (obj.y - s*real(phatfix)).^2 ./ (2*(obj.wvar + s^2*pvar));

                    % Combine to form output cost
                    ll = ls + 0.5*(real(Axhat - phatfix)).^2./pvar;
                end;

            else
                % Output cost is simply the log likelihood
                ll = -0.5*((obj.y-s*real(Axhat)).^2)./wvar1;
            end 
        end
        
        function S = numColumns(obj)
            %Return number of columns of Y
            S = size(obj.y,2);
        end
        
        % Generate random samples from p(y|z)
        function y = genRand(obj, z)
            y = sqrt(obj.wvar).*randn(size(z)) + obj.scale.*z;
        end
    end
    
end

