function sigmapow_RL = gammaRLfunc_AMP_ver2(rho,theta,mu_RL_prev,sigmapow_RL_prev,normpdfRL_matrix,theta_plus_sigmapow,q,sigmapow0,Z,mu_RL)


temp=theta.*mu_RL_prev+rho.*sigmapow_RL_prev;

sigmapow_RL...
=(  (1-q)*(  temp.^2                  + theta.* sigmapow_RL_prev              .* theta_plus_sigmapow(:,1))./ theta_plus_sigmapow(:,1).^2 .*normpdfRL_matrix(:,1)...
    +  q *( (temp + rho*sigmapow0).^2 + theta.*(sigmapow_RL_prev + sigmapow0) .* theta_plus_sigmapow(:,2))./ theta_plus_sigmapow(:,2).^2 .*normpdfRL_matrix(:,2)...
)./Z - (mu_RL.^2);
end