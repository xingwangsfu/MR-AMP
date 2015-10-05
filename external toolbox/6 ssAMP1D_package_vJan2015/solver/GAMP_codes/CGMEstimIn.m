classdef CGMEstimIn < EstimIn
    % CGMEstimIn:  Complex-valued Gaussian Mixture scalar input estimation function
    
    properties 
        omega; % Weights
        theta;  % Means 
        phi;   % Variances 
    end
    
    methods
        % Constructor
        % omega: weight or probability of a given model component (need not be non-normalized to 1)
        % theta: Mean of a given component
        % phi: Variance of a given component
        %
        % The arguments are 3-dimensional matrices, the first two dimensions correspond to
        % the dimensions of the matrix to be estimated. 
        % e.g. if the arguments are size (1000,5,3), the input vector is 1000x5 and its pdf has 3
        % Gaussian components
        function obj = CGMEstimIn(omega, theta, phi)
            obj = obj@EstimIn;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.omega = omega;
                obj.theta = theta;
                obj.phi = phi;            

                % normalize omega so sum(obj.omega, 3) ==1
                L = size(obj.omega,3);
                obj.omega = obj.omega ./ repmat(sum(obj.omega, 3), [1, 1, L]);
            end
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = sum(obj.theta.*obj.omega,3);
            var0  = sum(obj.omega .* (obj.phi + ...
                abs(obj.theta).^2), 3) - abs(mean0).^2;
            valInit = 0;
        end 


        function [Uhat, Uvar, NKL] = estim(obj, Rhat, Rvar)

            %Get the number of mixture components
            L = size(obj.omega,3);
            
            %Grab the signal dimension
            [N, T] = size(Rhat);

            %Expand scalar estimator if needed
            omega = obj.omega;
            theta = obj.theta;
            phi = obj.phi;
            if ((size(omega,1)==1)&(N>1))
              omega = repmat(omega,[N,1,1]);
              theta = repmat(theta,[N,1,1]);
              phi = repmat(phi,[N,1,1]);
            end
            if ((size(omega,2)==1)&(T>1))
              omega = repmat(omega,[1,T,1]);
              theta = repmat(theta,[1,T,1]);
              phi = repmat(phi,[1,T,1]);
            end

            %Preallocate storage
            gamma = zeros(N,T,L); alpha = zeros(N,T,L);
            beta = zeros(N,T,L); nu = zeros(N,T,L);

            for i = 1:L
               beta(:,:,i) = phi(:,:,i) + Rvar + eps;
               alpha(:,:,i) = abs(Rhat-theta(:,:,i)).^2./beta(:,:,i);
               gamma(:,:,i) = (Rhat.*phi(:,:,i) + theta(:,:,i).*Rvar)./beta(:,:,i);
               nu(:,:,i) = Rvar.*phi(:,:,i)./beta(:,:,i);
            end

            lik = zeros(N,T,L);
            for i = 1:L
                lik = lik + repmat(omega(:,:,i),[1 1 L])./omega...
                    .*beta./repmat(beta(:,:,i),[1 1 L])...
                    .*exp((alpha-repmat(alpha(:,:,i),[1 1 L])));
            end

            Uhat = sum(gamma./lik,3);
            Uvar = sum((nu + abs(gamma).^2)./lik,3) - abs(Uhat).^2;
            
            % Compute the negative KL divergence            
            if (nargout >= 3)                            
                zeta = sum(omega.*exp(-alpha)./(pi*beta),3);
                zeta(zeta == 0) = eps;
                NKL = log(zeta)+ log(pi*Rvar)+(Uvar + abs(Uhat - Rhat).^2)./(Rvar);
            end

        end
        
        % Generate random samples
        function x = genRand(obj, outSize) 

            if (size(obj.omega,1)~=1)||(size(obj.omega,2)~=1)
                error('genRand() implemented only for scalar CGMEstimIn');
            end 
            
            if isscalar(outSize)
                row = outSize;
                col = 1;
            else
                row = outSize(1);
                col = outSize(2);
            end

            L = size(obj.omega,3);
            omega = squeeze(obj.omega(1,1,:));
            theta = squeeze(obj.theta(1,1,:));
            phi = squeeze(obj.phi(1,1,:));

            dummy = [0;cumsum(omega)];
            dummy2 = rand(row,col);
            dummy3 = zeros(row,col,L);
            for i = 1:L
                dummy3(:,:,i) = ((dummy2>=dummy(i))&(dummy2<dummy(i+1)))...
                    .*(theta(i) + sqrt(phi(i)/2)...
                    .*(randn(row,col)+1i*randn(row,col)));
            end
            x = sum(dummy3,3);
        end
        
        % Computes the likelihood p(y) for y = x + v, v = N(0,Yvar)
        function py = plikey(obj,Y,Yvar)
            
            L = size(obj.omega,3);
            [M, T] = size(Y);
            lik = zeros(M,T,L);
            
            for i = 1:L
                lik(:,:,i) = obj.omega(:,:,i).*exp(-1./(obj.phi(:,:,i)+Yvar)...
                    .*abs(Y-obj.theta(:,:,i)).^2)./(pi*(obj.phi(:,:,i)+Yvar));
            end
            
            lik(isnan(lik)) = 0.999;
            py = sum(lik,3);
        end
        
        % Computes the log-likelihood, log p(Y(i,j)), for Y = X + V, where 
        % p(X(i,j)) = sum_k omega(i,j,k)*CN(theta(i,j,k), phi(i,j,k)) and 
        % p(V(i,j)) = CN(0, Yvar(i,j))
        function logpy = loglikey(obj, Y, Yvar)
            logpy = log(obj.plikey(Y, Yvar));
        end

    end
    
end

