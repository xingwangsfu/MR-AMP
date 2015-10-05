function [X_vec] = patch_unite(Y,Params)

M = Params.M;
N = Params.N;
n = Params.n;
m = Params.m;
ids_vec = Params.ids_vec;
blocksize = Params.blocksize;
stepsize = Params.stepsize;
blocknum = Params.blocknum;
M_dim = Params.M_dim;
N_dim = Params.N_dim;

X = zeros(M,N);

for i = 1:blocknum
    ids{1} = [ids_vec{i}(1):ids_vec{i}(1)+m-1];
    ids{2} = [ids_vec{i}(2):ids_vec{i}(2)+n-1];
    tmp = reshape(Y(:,i),m,n);
    X(ids{1},ids{2}) = X(ids{1},ids{2})+tmp;
end

cnt = countcover(size(X),blocksize,stepsize);
X = X./cnt;
zero_pos = find(cnt==0);
X_vec = X(:);
X_vec(zero_pos) = [];
X_vec = reshape(X_vec,M_dim,N_dim);