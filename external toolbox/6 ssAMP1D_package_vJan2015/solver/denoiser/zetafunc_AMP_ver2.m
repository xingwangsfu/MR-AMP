function Z=zetafunc_AMP_ver2(normpdf_matrix,C,q)

    Z= (1-q)^2  .* C(:,1) .*  normpdf_matrix(:,1) +...
       q*(1-q)  .* C(:,2) .*  normpdf_matrix(:,2)+...
       q*(1-q)  .* C(:,3) .*  normpdf_matrix(:,3)+...
       q^2      .* C(:,4) .*  normpdf_matrix(:,4);
   
%    TH_Z=1e-7;
%        Z(find(Z<TH_Z))=TH_Z;  %To prevent NaN   
%  
end
