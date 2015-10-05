classdef SoftThreshEstimIn < EstimIn
    % SoftThreshEstimIn:  Inputs a soft thresholding scalar input function.
    % Allows GAMP to be used to solve "min_x 1/2/var*norm(y-A*x,2)^2 + lambda*norm(x,1)". 
    
    properties
        %lambda = the gain on the ell1 term in the MAP cost. 
	%The soft threshold is set according to the expression thresh = lambda * mur;
        lambda;
        maxSumVal = true;   % Max-sum GAMP (true) or sum-product GAMP (false)?
    end
    
    methods
        % Constructor
        function obj = SoftThreshEstimIn(lambda, maxSumVal)
            obj = obj@EstimIn;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.lambda = lambda;
                if nargin >= 2 && ~isempty(maxSumVal) && isscalar(maxSumVal)
                    obj.maxSumVal = logical(maxSumVal);
                end
            end
        end
        
        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(~)
	    %Set these to reasonable constants
            mean0 = 0; 
            var0  = 5e-4;
            valInit = -inf;
        end
        
        % Carry out soft thresholding
        function [xhat,xvar,val] = estim(obj, rhat, rvar)
            if ~obj.maxSumVal
                % Compute sum-product GAMP updates
                
                % To avoid numerical problems (0/0) when evaluating
                % ratios of Gaussian CDFs, impose a firm cap on the
                % maximum value of entries of rvar
                rvar = min(rvar, 700);
                
                % *********************************************************
                % Begin by computing various constants on which the
                % posterior mean and variance depend
                sig = sqrt(rvar);                       	% Gaussian prod std dev
                muL = rhat + obj.lambda.*rvar;          	% Lower integral mean
                muU = rhat - obj.lambda.*rvar;          	% Upper integral mean
                muL_over_sig = muL ./ sig;
                muU_over_sig = muU ./ sig;
                cdfL = normcdf(-muL_over_sig);              % Lower cdf
                cdfU = normcdf(muU_over_sig);               % Upper cdf
                cdfRatio = cdfL ./ cdfU;                    % Ratio of lower-to-upper CDFs
                SpecialConstant = exp( (muL.^2 - muU.^2) ./ (2*rvar) ) .* ...
                    cdfRatio;
                NaN_Idx = isnan(SpecialConstant);        	% Indices of trouble constants
                
                % For the "trouble" constants (those with muL's and muU's
                % that are too large to give accurate numerical answers),
                % we will effectively peg the special constant to be Inf or
                % 0 based on whether muL dominates muU or vice-versa
                SpecialConstant(NaN_Idx & (-muL >= muU)) = Inf;
                SpecialConstant(NaN_Idx & (-muL < muU)) = 0;
                
                % Compute the ratio normpdf(a)/normcdf(a) for
                % appropriate upper- and lower-integral constants, a
                RatioL = 2/sqrt(2*pi) ./ erfcx(muL_over_sig / sqrt(2));
                RatioU = 2/sqrt(2*pi) ./ erfcx(-muU_over_sig / sqrt(2));
                
                % Now compute the first posterior moment...
                xhat = (1 ./ (1 + SpecialConstant.^(-1))) .* ...
                    (muL - sig.*RatioL) + (1 ./ (1 + SpecialConstant)) .* ...
                    (muU + sig.*RatioU);
                
                % ...and second central posterior moment
                varL = rvar .* (1 - RatioL.*(RatioL - muL_over_sig));
                varU = rvar .* (1 - RatioU.*(RatioU + muU_over_sig));
                meanL = muL - sig.*RatioL;
                meanU = muU + sig.*RatioU;
                SecondMoment = (1 ./ (1 + SpecialConstant.^(-1))) .* ...
                    (varL + meanL.^2) + (1 ./ (1 + SpecialConstant)) .* ...
                    (varU + meanU.^2);
                xvar = SecondMoment - xhat.^2;
                % *********************************************************
                
                % Lastly, compute negative KL divergence:
                % \int_x p(x|r) log(p(x)/p(x|r)) dx
                
%                 % Old way of computing.  It handles difficult cases
%                 % incorrectly.
%                 NormConL = obj.lambda/2 .* ...              % Mass of lower integral
%                     exp( (muL.^2 - rhat.^2) ./ (2*rvar) ) .* cdfL;
%                 NormConU = obj.lambda/2 .* ...              % Mass of upper integral
%                     exp( (muU.^2 - rhat.^2) ./ (2*rvar) ) .* cdfU;
%                 NormCon = NormConL + NormConU;      % Posterior normaliz. constant recip.
%                 NormCon(isnan(NormCon)) = 1;        % Not much we can do for trouble ones
%                 NormCon(NormCon == Inf) = 1;
                
                if (nargout >= 3)  
                    %Calculate lower and Upper non-shared integration factors
                    CL = erfcx(muL_over_sig/sqrt(2));
                    CU = erfcx(-muU_over_sig/sqrt(2));
                    % The individual terms can still be infinite. For these
                    % cases use the approximation erfc(x) = 1/sqrt(pi)/x for
                    % large x
                    I = find(isinf(CL));
                    CL(I) = 1./(sqrt(pi/2)*abs(muL_over_sig(I)));
                    I = find(isinf(CU));
                    CU(I) = 1./(sqrt(pi/2)*abs(muU_over_sig(I)));

                    % Calculate the log scale factor
                    logNormCon = log(obj.lambda/4) - rhat.^2./(2*rvar) + log(CL + CU);

                    val = logNormCon + ...
                            0.5*(log(2*pi*rvar) + ...  
                            (xvar + (xhat - rhat).^2)./rvar);
                    
                    %Old fix, no longer needed
                    %val(val == -Inf) = -1e4;    % log(NormCon) can == -Inf
                end
            else
                % Compute max-sum GAMP updates
                
                %Compute the thresh
                thresh = obj.lambda .* rvar;
                
                %Estimate the signal
                xhat = max(0,abs(rhat)-thresh) .* sign(rhat);
                
                %Estimate the variance
                %xvar = rvar .* (mean(double(abs(xhat) > 0))*ones(size(xhat)));
                xvar = rvar .* (abs(xhat) > 0);
                
                if (nargout >= 3)  
                    %Output negative cost
                    %val = -1*obj.lambda*abs(rhat);
                    val = -1*obj.lambda.*abs(xhat);	% seems to work better
                end
            end
        end
        
        % Computes p(y) for y = x + v, with x ~ p(x), v ~ N(0,yvar)
        function py = plikey(obj, y, yvar)
            mu = y;
            sig2 = yvar;
            sig = sqrt(sig2);                           % Gaussian prod std dev
            muL = mu + obj.lambda.*sig2;                % Lower integral mean
            muU = mu - obj.lambda.*sig2;                % Upper integral mean
            muL_over_sig = muL ./ sig;
            muU_over_sig = muU ./ sig;
            erfsum = erfcx((1/sqrt(2))*muL_over_sig) + ...
                erfcx((-1/sqrt(2))*muU_over_sig);
            erfsum(erfsum == Inf) = realmax;
            py = obj.lambda/4 .* exp(-(y.^2) ./ (2*yvar)) .* erfsum;
        end
        
    end
    
end

