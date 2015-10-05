function [tau_new] = appro_gd(delta_N, kappa, alpha, beta, mu_tilde, m, ell_zero, Flag, N, r_tilde)

index = 0;
while Flag ==1 && index <= 10 
    Flag = 0;
    tau_new = 0;
    for i = 1:m
        tau_old = tau_new;
        df_r_tilde = (abs(r_tilde(tau_old+delta_N)) - abs(r_tilde(tau_old)))/delta_N;
        delta_tau = -1*df_r_tilde;
        ell = ell_zero;
        while abs(r_tilde(tau_old+ell*delta_tau)) > abs(r_tilde(tau_old)) + alpha*ell*delta_tau*df_r_tilde
            ell = beta*ell;
        end
        tau_new = (tau_old + ell*delta_tau);
        if abs(abs(r_tilde(tau_new))/N-mu_tilde)/mu_tilde < kappa ||   tau_new <0
            ell_zero = ell_zero/2; Flag = 1; index = index+1;
            break;
        end
    end
    if i == m 
        Flag = 0;
        tau_new = abs(tau_new);
    end
end
tau_new = abs(tau_new);
            