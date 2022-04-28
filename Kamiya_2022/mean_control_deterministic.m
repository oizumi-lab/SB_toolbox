function [ret, M_T0, Phi_T0] = mean_control_deterministic(A, BtB, mu_0, mu_T, t0, tf, delta_t)
n = size(A,1);
stp_n = fix( (tf-t0)/delta_t ) + 1;



M_t0 = zeros(n,n,stp_n);
M_t_t0 = zeros(n,n,stp_n-1);
for ii = 2:stp_n
    M_t_t0(:,:,ii-1) = A * M_t0(:,:,ii-1) + M_t0(:,:,ii-1)*A' + BtB;
    M_t0(:,:,ii) = M_t0(:,:,ii-1) + M_t_t0(:,:,ii-1)*delta_t;

end
M_T0 = M_t0(:,:,end);    
    
Phi_Ts_mat = Phi_Ts(A, t0, tf, delta_t);
Phi_t0_mat = zeros(n,n,stp_n);
for kk = 1:stp_n
    Phi_t0_mat(:,:,kk) = inv(Phi_Ts_mat(:,:,kk))*Phi_Ts_mat(:,:,1);
end
Phi_T0 = Phi_t0_mat(:,:,end);

m_tmp = (mu_T - Phi_T0 * mu_0);
ret = m_tmp' * inv(M_T0) * m_tmp;
end