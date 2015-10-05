% this file is to test the phase transition curve of TVAMP for 1-D
% piecewise-constan signal
clear all;
N = 628;
d=2;
delta_cand = linspace(0.1,0.4,10);
% delta_cand = 0.5;
% rho_cand = 0.5/d;
stop_tol = (1e-12)^2;
Up_matrix = zeros(N,N/d);
%downFactor = N/n;
if mod(d,1) ==0
    for i = 1:size(Up_matrix,2)
        Up_matrix((i-1)*d+1:i*d,i) = 1;
    end
end
maxiter = 2000;
epsilon = 0.05;
SNR = 60; % 60dB
for deltaIdx = 1:length(delta_cand)
    delta = delta_cand(deltaIdx);
    M = floor(delta*N);
    K = floor(epsilon*(N));
    count = 0;
    count_LR = 0;
    parfor index = 1:100
        x_LR=Gaussian_1DFD_gen(N/d,d*epsilon,1);
        x_mtx = repmat(x_LR',d,1);
        x = x_mtx(:);
        Psi = randn(M,N)/sqrt(N);
        normalized_vector = sqrt(sum(Psi.*Psi,1));
        normalized_matrix = repmat(normalized_vector,M,1);
        A = Psi./normalized_matrix;
        AT = A';
        y = A*x;
        sigma = sqrt(norm(y)^2/(M*10^(SNR/10)));
        y = y+sigma*randn(M,1);
        
        % HR reconstruction
        %
        [x_rec,tau]=solve_TVAMP_mat_ver4(A,y,maxiter,x,stop_tol);
        NMSE(deltaIdx,index) = norm(x_rec-x)^2/norm(x)^2;
        
        % equivalent measurement matrix for LR-AMP
        A_LR = A*Up_matrix/sqrt(d);
        [x_rec_LR_un,tau] = solve_TVAMP_mat_ver4(A_LR,y,maxiter,sqrt(d)*x_LR,d*stop_tol);
        NMSE_LR(deltaIdx,index) = norm(x_rec_LR_un/sqrt(d)-x_LR)^2/norm(x_LR)^2;
    end
    deltaIdx
end
