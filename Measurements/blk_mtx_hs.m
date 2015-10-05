function y = blk_mtx_hs(z,Psi,N,S)


M = size(Psi,1);
z = reshape(z,[S N N]);
y_mtx = zeros(S,M);
for i = 1:S
    alpha_tmp = z(i,:,:);
    y_mtx(i,:) = Psi*alpha_tmp(:);
end

y = y_mtx(:);


