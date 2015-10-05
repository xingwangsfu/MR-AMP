% This "inverts" the function Axhat = obj.estim(phat;pvar) in that
% it solves for the input phat that leads to a specified output Axhat

function [phat,zhat,zvar] = estimInvert(obj,Axhat,pvar,opt)

    if (nargin<4)||(~isfield(opt,'debug')), opt.debug = false; end
    if (nargin<4)||(~isfield(opt,'phat0')), opt.phat0 = Axhat; end
    if (nargin<4)||(~isfield(opt,'alg')), opt.alg = 1; end
    if (nargin<4)||(~isfield(opt,'maxIter')), opt.maxIter = 50; end;
    if (nargin<4)||(~isfield(opt,'tol')), opt.tol = 1e-4; end;
    if (nargin<4)||(~isfield(opt,'stepsize')), opt.stepsize = 0.5; end;
    if (nargin<4)||(~isfield(opt,'regularization')), opt.regularization = 0; end;

    if opt.debug, 
      phat_ = nan(length(Axhat(:)),opt.maxIter); 
      zhat_ = nan(length(Axhat(:)),opt.maxIter); 
    end;

    % Test quality of initialization
    phat = opt.phat0;
    [zhat,zvar] = obj.estim(phat,pvar);
    residual = Axhat-zhat;
    NR0 = norm(residual);
    zhat0 = zhat;

    % Iteratively improve fixed-point estimate
    r=0;
    NR = NR0;
    gamma = opt.stepsize;
    phatOld = inf(size(phat));
    while (r<opt.maxIter)&&(NR>opt.tol*norm(zhat))
    %while (r<opt.maxIter)&&(norm(phat-phatOld)>opt.tol*norm(phat))

	r = r+1;
	if opt.debug
            phat_(:,r) = phat(:);
            zhat_(:,r) = zhat(:);
	end

        % Update estimate
	phatOld = phat;
	switch opt.alg
	    case 0 % simple method
                phat = phatOld + gamma*residual; 
	    case 1 % approximate newton's method
	        gradient = zvar./pvar;
                phat = phatOld + gamma*residual.*gradient./(gradient.^2 + ...
		 	+ opt.regularization); 
	end
	 
        % Test quality 
        [zhat,zvar] = obj.estim(phat,pvar);
        residual = Axhat-zhat;
	NR = norm(residual);
    end

    % Check if progress was made
    if (NR0<0.9*NR)
       warning('No progress made: Adjust estimInvert options in *EstimOut.logScale') 
    end

    % Present debug information
    if opt.debug
        % report iterations and error
	fprintf(1,'iterations=%3g\n',r);
	fprintf(1,'norm(zhat-Axhat)=%5.5g\n',...
		norm(zhat-Axhat))
	fprintf(1,'norm(zhat-Axhat)/norm(zhat)=%5.5g\n',...
		norm(zhat-Axhat)/norm(zhat))
	fprintf(1,'norm(phat-phatOld)/norm(phat)=%5.5g\n\n',...
		norm(phat-phatOld)/norm(phat))

        % plot example trajectories
        figure(100)
	num_coef = min(100,size(Axhat,1));
	mm = 1:num_coef;
	subplot(211)
            plot([1:r],real(phat_(mm,1:r))')
	    grid on
	    ylabel('real(phat)')
	    xlabel('iteration')
	title(['First ',num2str(num_coef),' estimInvert coefficient trajectories'])
	subplot(212)
	    Axhat_ = Axhat(:);
	    semilogy([1:r],abs(zhat_(mm,1:r)-Axhat_(mm)*ones(1,r))')
	    grid on
	    ylabel('abs(zhat-Axhat)')
	    xlabel('iteration')
	drawnow
	pause
    end

