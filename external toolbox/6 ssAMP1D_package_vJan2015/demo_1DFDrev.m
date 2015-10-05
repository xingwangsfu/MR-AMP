%------------------------------------------------------------------------%
% Filename: demo_1DFDrev.m
% This file is a testbench for an exeamplary comparison among 
% algorithms for the piecewise-constant recovery problem.  
% 
% The list of piecewise-constant recovery algorithm in this comparison.
% [0] ssAMP-1D: the corresponding paper -  Jaewook Kang, Heung-No Lee, and Kiseon Kim,  
%     "One-dimensional  Piecewise-Constant Signal Recovery via  Spike-and-Slab Approximate Message-Passing,"  
%      to appear in proc. of the 48th Asilomar Conference (Pacific Grove, CA), Nov. 2014.
% [1] EFLA : the corresonding paper - J. Liu, L. Yuan, and J. Ye. 
%    ¡°An efficient algorithm for a class of fused lasso problems,¡± 
%     proc of ACM SIGKDD Conference on Knowledge Discovery and Data Mining, 
%     2010.
% [2] TV-AMP: the corresonding paper - D. L. Donoho, I. Johnstone, and A. Montanari, 
%     ¡°Accurate prediction of phase transitions in compressed sensing via 
%      a connection to minimax denoising, ¡± IEEE Trans. Inform. Theory, 
%      vol. 59, no. 6, pp. 3396-3433, June 2013.
% [4] Chambolle-Pock: the corresonding paper - A. Chambolle, T. Pock, 
%    ¡°A first-order primal-dual algorithm for convex
%      problems with applications to imaging,¡± J. Math. Imag. Vis., vol. 40, pp.
%      120-145, May 2011
% [5] GrAMPA: the corresonding paper - M. Borgerding and P. Schniter, 
%     ¡°Generalized approximate message passing for the cosparse analysis
%      model,¡± avabilable at ArXiv:1312.3968v1 [cs.IT], Dec. 2013
% written by Jaewook Kang
% Final updated at Jan, 2015
%------------------------------------------------------------------------%
clc
clear all
close all

%Handle random seed

if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end;

if 1
    savedState = defaultStream.State;
    save random_state.mat savedState;
else
    load random_state.mat
end

defaultStream.State = savedState;

% put key subdirectories in path if not already there
path(path, './solver/flsa');
path(path, './solver/denoiser');
path(path, './solver/GAMP_codes');
path(path, './solver');
path(path, './etc_tool');
path(path,genpath(pwd));

disp('%----------------------------------------------------------------------------------------------------------------------%');
disp('% Spike-and-Slab approximate message-passing recovery for one-dimensional piecewise-constant signals');
disp('% Copyright@Jaewook Kang with CSNL lab in GIST, Republic of Korea.');
disp('%');
disp('% Written by Jaewook Kang, Phd student in GIST-DIC, jwkkang@gist.ac.kr');
disp('% final update 2014 Jan');
disp('%');
disp('% This work was supported by Do-Yak  (NO.2013-035295),and Leading Foreign Research Institute (MT-IT)  (2009-00422) Progrem');
disp('% Special thank to prof. Schniter Phil (Ohio state University) for help to construct the experiment');
%--------------------- Problem dimension setting -------------------------%
alpha=0.5;% undersampling ratio M/N %0.5 / 0.05 / 0.8
KM=0.05;
N=1000; % the signal vector dimension %25^2 
M=round (N*alpha); % the measurement vector dimension
%  Additive noise variance
Delta=1e-5;
%---------------------- Choose solvers --------------------%
usessAMP1D = true;
useEFLA = true;
useTVAMP = true;
useCP = true;
useGrAMPAsnipe = true;
useGrAMPAbg = true;
%---------------------- Measurement matrix generation --------------------%
H=randn(M,N)/sqrt(M);% 1) Gaussian matrices        
%H=(2*(rand(M,N)<0.5)-1)/sqrt(M);% 2) Bernoulli matrices h \in {1,-1}
%H=randint2(M,N,[-1,1])/sqrt(dilution*M);% 3) Ternary matrices h \in {0,1,-1}


%----------------Piecewise constant signal generation --------------------%
% These parameters are applied to the spike-and-slab potentional function
% in the ssAMP1D algorithm
q      = KM*alpha;% signal sparsity
sigma0 = 1;% variance for difference between neigboring values


X=Gaussian_1DFD_gen(N,q,sigma0);% Gaussian signal: X_i - X_i-1 \in Gaussian
%X=Ternary_1DFD_gen(N,q,sigma0);% Ternary signal: X_i - X_i-1 \in {0,1,-1}
%clf; plot(X); figure(gcf); return

disp('%----------------------------------------------------------------------------------------------------------------------%');
disp('<Experiment condition>')
disp(sprintf('Signal length, N = %d',N));
disp(sprintf('Piecewise-constancy, K/M = %8.3f',q/alpha));
disp(sprintf('Undersampling ratio, M/N  = %8.3f',M/N));
disp(sprintf('Noise variance, delta  = %8.3f',Delta));
disp('%----------------------------------------------------------------------------------------------------------------------%');
disp('<The number of iteration expended>')
%--------- Measurement generation by random linear projection-------------%

noise=sqrt(Delta)*randn(M,1);
y=H*X+noise;
% Common parameters
maxiter = 1000;
wvar = max(Delta,1e-10);
SNR = (mean(abs(y).^2)-wvar)/wvar; 
iteration_stopping_tol = (max(1e-12,min(1e-4,(1/SNR))))^2;  
%------------parameters for CP--------------------%
CP_lambda=0.5;
%--------------------- option setting for EFLA solver --------------------%
rho=0;          % the regularization parameter;it is a ratio between (0,1), if .rFlag=1
opts_EFLA=[];
opts_EFLA.init=2;        % starting from a zero point
% termination criterion
% opts_EFLA.tFlag=5;       % run .maxIter iterations
opts_EFLA.tFlag=3;       % iteration termination when norm( x_i - x_{i-1}, 2) <= .tol


opts_EFLA.maxIter=maxiter;   % maximum number of iterations
opts_EFLA.nFlag=0;       % without normalization
% regularization
opts_EFLA.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
opts_EFLA.fusedPenalty=0.1;
opts_EFLA.lFlag=0; % line search
%----------------------- option setting for TVAL3 ------------------------%
opts_TVAL3.mu =  2^7;
opts_TVAL3.beta = 2^4;
opts_TVAL3.mu0 = opts_TVAL3.mu;       % trigger continuation shceme
opts_TVAL3.beta0 = opts_TVAL3.beta;    % trigger continuation shceme
opts_TVAL3.maxit=maxiter; 
opts_TVAL3.maxcnt =maxiter;
opts_TVAL3.tol_inn =  sqrt(iteration_stopping_tol);
opts_TVAL3.tol =1e-6;
opts_TVAL3.init=0; % start from the zero point
opts_TVAL3.TVnorm = 1; % Anisortopic setting
opts_TVAL3.TVL2=true;  % The TV/L2+ model
%------------------------- Configuring GrAMPA ----------------------------%
% 1D TV dictionary 
I = speye(N-1);
P = spalloc(N-1,1,0);
scaling = sqrt(norm(H,'fro')^2/M); % to equalize row norms of H and D
D = scaling/sqrt(2) *( [P I] - [I P] ); % note sparse and properly scaled
% linear transform exploiting nonsparsity of H and sparsity of D
% GrampaLinTrans = LinTransConcat({MatrixLinTrans(H);MatrixLinTrans(D)});
GrampaLinTrans = MatrixLinTrans([H;D]);

% options that seem to work well for 1D TV dictionaries
 GrampaOptions = GampOpt;
 GrampaOptions.legacyOut = false;
 GrampaOptions.xvar0 = (norm(y)^2-Delta*M)/norm(H,'fro')^2; % guess at signal variance
 xvar0big = 100*GrampaOptions.xvar0; % much larger than signal variance 
 
 GrampaOptions.tol = sqrt(iteration_stopping_tol); 
 GrampaOptions.nit = maxiter; % maximum number of iterations
 GrampaOptions.varNorm = false; % turn off internal normalization
 GrampaOptions.zvarToPvarMax = inf; % do not clip 
 GrampaOptions.uniformVariance = true; % off since applied externally
 GrampaOptions.xvarMin = 0; % no minimum rvar
 GrampaOptions.pvarMin = 0; % no minimum pvar
%  GrampaOptions.pvarStep = true; % apply damping to pvar 
  GrampaOptions.pvarStep = false; % apply damping to pvar 
 GrampaOptions.adaptStep = false; % leave off 
 %-----------------------------------------------------------
%  GrampaOptions.adaptStepBethe = true; % use Bethe version 
%  GrampaOptions.stepWindow = 10; % adaptive stepsize window 
%  GrampaOptions.stepIncr = 1.1; % stepsize increase rate
%  GrampaOptions.stepDecr = 0.5; % stepsize decrease rate
 %--------------------------------------------------------
 GrampaOptions.stepMax = 1.0; % maximum stepsize: 1.0 for speed, 0.5 for robustness
 GrampaOptions.stepMin = 0.25; % minimum stepsize
 GrampaOptions.step = GrampaOptions.stepMax; % initial stepsize

% trivial prior
% GrampaEstimIn = NullEstimIn(0,1);
GrampaEstimIn = AwgnEstimIn(0,xvar0big); % very mild regularization, for stability

% AWGN measurements 
MeasEstimOut = AwgnEstimOut(y,Delta+eps);

% SNIPE regularization 
omega = log((1-q)/q);  % SNIPE regularization _roughly_ tuned to q
AnaEstimOut_snipe = SNIPEstim(omega); % SNIPE denoising
GrampaEstimOut_snipe = EstimOutConcat({MeasEstimOut;AnaEstimOut_snipe},[M,N-1]);

% spike-and-slab regularization 
AnaEstimOut_ss = AwbgnEstimOut(0,(scaling^2/2)*sigma0^2,q); % spike-and-slab
GrampaEstimOut_ss = EstimOutConcat({MeasEstimOut;AnaEstimOut_ss},[M,N-1]);
%--------------------------- Signal recovery ----------------------------%

normX_sqr=norm(X)^2;

% ssAMP1D solving [0] (Proposed)
if usessAMP1D
tstart=tic;
[estX_fusedAMP,varX,end_iter,theta]= solve_ssAMP_mat_ver5(H,y,3*sigma0,q,Delta,maxiter,iteration_stopping_tol);
telapsed_fusedMP=toc(tstart);
MSE_fusedMP= norm(estX_fusedAMP-X)^2/normX_sqr;
end


% EFLA solving [1]
if useEFLA
% EFLA considers the stopping criterion with " norm( x_i - x_{i-1}, 2) <=.tol" 
% therfore,
opts_EFLA.tol=sqrt(iteration_stopping_tol*normX_sqr); 
    
tstart=tic;
[estX_fusedlasso, funVal1, ValueL1]= fusedLeastR(H, y, rho, opts_EFLA);
telapsed_fusedlasso=toc(tstart);
MSE_fusedlasso= norm(estX_fusedlasso-X)^2/normX_sqr;
end

% TV-AMP solving [2]
if useTVAMP
tstart=tic;
[estX_TVAMP,tau]=solve_TVAMP_mat_ver4(H,y,maxiter,X,iteration_stopping_tol);
telapsed_TVAMP=toc(tstart);
MSE_TVAMP= norm(estX_TVAMP-X)^2/normX_sqr;
end


% chambolle-pock TV [4]
if useCP
% x_init=H'*y;
x_init=zeros(N,1);
 L=max(svd(H));
tstart=tic;
estX_CP= solve_chambolle_pock_TV_ver2(H,y,x_init, CP_lambda,iteration_stopping_tol,maxiter,L);
telapsed_CP=toc(tstart);
MSE_CP=norm(estX_CP-X)^2/normX_sqr;
end


% GrAMPA estimate with SNIPE regularization [5]
if useGrAMPAsnipe
tstart=tic;
estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut_snipe,GrampaLinTrans,GrampaOptions);
telapsed_GrAMPAsnipe=toc(tstart);
estX_GrAMPAsnipe=estFin1.xhat;
MSE_GrAMPAsnipe=norm(estX_GrAMPAsnipe-X)^2/normX_sqr;
disp(sprintf('Number of GrAMPAsnipe iterations = %d',estFin1.nit));
end


% GrAMPA estimate with spike-and-slab regularization [5]
if useGrAMPAbg
tstart=tic;
estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions);
telapsed_GrAMPAbg=toc(tstart);
estX_GrAMPAbg=estFin1.xhat;
MSE_GrAMPAbg=norm(estX_GrAMPAbg-X)^2/normX_sqr;
disp(sprintf('Number of GrAMPAbg iterations = %d',estFin1.nit));
end

%--------------------------- Display ------------------------------------%
disp('%------------------------------------------------------------------------------------------%');
disp('<Recovery result>')
if usessAMP1D, disp(sprintf('ssAMP-1D: Nomalized MSE  = %8.7f',MSE_fusedMP)); end
if useEFLA, disp(sprintf('EFLA: Nomalized MSE  = %8.7f',MSE_fusedlasso)); end
if useTVAMP, disp(sprintf('TVAMP-FLSA: Nomalized MSE   = %8.7f',MSE_TVAMP)); end
if useCP, disp(sprintf('Chambolle-Pock L2TV:Nomalized MSE   = %8.7f',MSE_CP)); end
if useGrAMPAsnipe, disp(sprintf('GrAMPAsnipe: Nomalized MSE   = %8.7f',MSE_GrAMPAsnipe)); end
if useGrAMPAbg, disp(sprintf('GrAMPAbg: Nomalized MSE   = %8.7f',MSE_GrAMPAbg)); end
disp('%------------------------------------------------------------------------------------------%');
if usessAMP1D, disp(sprintf('Running Time of ssAMP-1D  = %8.7f sec',telapsed_fusedMP)); end
if useEFLA, disp(sprintf('Running Time of EFLA = %8.7f sec',telapsed_fusedlasso)); end
if useTVAMP, disp(sprintf('Running Time of TVAMP-FLSA= %8.7f sec',telapsed_TVAMP)); end
if useCP, disp(sprintf('Running Time of Chambolle-Pock L2TV= %8.7f sec',telapsed_CP)); end
if useGrAMPAsnipe, disp(sprintf('Running Time of GrAMPAsnipe= %8.7f sec',telapsed_GrAMPAsnipe)); end
if useGrAMPAbg, disp(sprintf('Running Time of GrAMPAbg= %8.7f sec',telapsed_GrAMPAbg)); end
disp('%------------------------------------------------------------------------------------------%');


figure(1); clf;
if usessAMP1D, subplot(3,2,1); plot(1:N,X,1:N,estX_fusedAMP,'rO');axis([1 N -max(abs(X))-1 max(abs(X))+1]); title('(a) ssAMP-1D recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useEFLA, subplot(3,2,2); plot(1:N,X,1:N,estX_fusedlasso,'rO');axis([1 N -max(abs(X))-1 max(abs(X))+1]); title('(b) EFLA recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useTVAMP, subplot(3,2,3); plot(1:N,X,1:N,estX_TVAMP,'rO');axis([1 N -max(abs(X))-1 max(abs(X))+1]); title('(c) TV-AMP recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useCP, subplot(3,2,4); plot(1:N,X,1:N,estX_CP,'rO');axis([1 N -max(abs(X))-1 max(abs(X))+1]); title('(e) Chambolle-Pock recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useGrAMPAsnipe, subplot(3,2,5); plot(1:N,X,1:N,estX_GrAMPAsnipe,'rO');axis([1 N -max(abs(X))-1 max(abs(X))+1]); title('(f) GrAMPAsnipe recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useGrAMPAbg, subplot(3,2,6); plot(1:N,X,1:N,estX_GrAMPAbg,'rO');axis([1 N -max(abs(X))-1 max(abs(X))+1]); title('(f) GrAMPAbg recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
box on;


% 
