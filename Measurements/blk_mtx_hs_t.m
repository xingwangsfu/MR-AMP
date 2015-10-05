function alpha = blk_mtx_hs_t(y,Psi,N,S)


M = size(Psi,1);
y = reshape(y,S,M);

alpha_mtx = zeros(S,N,N);
for i = 1:S
    tmp = Psi'*y(i,:)';
    tmp_mtx = reshape(tmp,N,N);
    alpha_mtx(i,:,:) = tmp_mtx;
end

alpha = alpha_mtx(:);