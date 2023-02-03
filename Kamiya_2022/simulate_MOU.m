function data_MOU = simulate_MOU(x0, A, C, TR, dt, T)

dim = size(A,1);
OU_mat = zeros(T/dt+1, dim);
OU_mat(1,:) = x0';

for ii = 2:T/dt+1
    OU_mat(ii,:) = ( OU_mat(ii-1,:)' + A * OU_mat(ii-1,:)' * dt +...
        sqrt(dt) * (C)^(1/2) * randn(dim,1) )';
end

data_MOU = downsample(OU_mat,TR/dt);
data_MOU = real(data_MOU);

end
