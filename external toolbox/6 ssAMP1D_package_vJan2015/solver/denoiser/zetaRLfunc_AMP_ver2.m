function Z = zetaRLfunc_AMP_ver2(normpdfRL_matrix,q)

    Z= (1-q).* normpdfRL_matrix(:,1) +   q .* normpdfRL_matrix(:,2);
% 
%    TH_Z=1e-7;
%        Z(find(Z<TH_Z))=TH_Z;  %To prevent NaN   
  
    
end