% this file is to test the empirical phase transition curve of
% multi-resolution compressed sensing

% generate the Haar wavelet transform
clear all;
L = 2; % level
%H =HaarTrans(L); % forward Haar transformation matrix

N = 1024;
n = 2^(L);

Psi = eye(N);
Psi = Psi';
factor = 2;
DeltaCand = 0.2;
RhoCand = 0.2;
%SuccessNum_LR = ones(30,30)*1e3;
% factor = 4;
for DeltaIndex = 1:length(DeltaCand)
    for RhoIndex = 1:length(RhoCand)
        delta = DeltaCand(DeltaIndex);
        rho = RhoCand(RhoIndex);
        epsilon = delta*rho;
        M = floor(delta*N);
        K = floor(epsilon*N);
        epsilon_LR = factor*epsilon;
        count = 0;
        count_LR = 0;
        if (factor == 4) % compute the thresholding parameter of soft-thresholding function using "minimax" rule
            [alpha,lambda] = optimal_threshold_DMM(delta);
            [alpha,lambda_LR] = optimal_threshold_DMM(factor^2*delta);
        else
            lambda = 0;
            lambda_LR = 0;
        end
        
        for MCindex = 1:100
            Phi = 1/sqrt(M)*randn(M,N);
            normalized_vector = sqrt(sum(Phi.*Phi,1));
            normalized_matrix = repmat(normalized_vector,M,1);
            Phi = Phi./normalized_matrix;
            x = zeros(N,1);
            chosenIdx = randperm(N); % index of sub-blocks
            x(chosenIdx(1:K))=randn(K,1);
            
            % measurement vector
            A = Phi*Psi;
            y = A*x;
            y = y + randn(size(y));
            % do AMP reconstruction
            x_LR_ref = x(1:N/factor);
            A_LR = A(:,1:N/factor);
            
            %    A_LR = Phi*Up_matrix/sqrt(factor);\
            Params.x_SI = 0; % no side information
            Params.T = 60; % number of iterations
            Params.tol = 1e-8; % tolerance
            Params.sigma_SI = 0;
            Params.mode = 'Auto'; % two modes here: auto and minimax, auto for parameterless AMP, minimax for DMM
            Params.SIorNot = 'Not'; % two modes: SI exists or not
            Params.x = zeros(N,1);
            [x_rec, pseudo_data, xi] = SI_AMP(y, A, A', N, lambda, Params);
            
            [alpha_L2H] = STAMP_L2H(y,A,A',N,lambda,pseudo_data,xi);
            %             [x_rec_LR] = SI_AMP(y, A_LR, A_LR', N/factor, lambda_LR);
            %             if norm(x_rec-x)^2/norm(x)^2 < 1e-6
            %                 count = count + 1;
            %             end
            %
            %             if norm(x_rec_LR-x_LR_ref)^2/norm(x_LR_ref)^2 < 1e-6
            %                 count_LR = count_LR + 1;
            %             end
        end
        
        SuccessNum(DeltaIndex,RhoIndex) = count;
        SuccessNum_LR(DeltaIndex,RhoIndex) = count_LR;
    end
end
