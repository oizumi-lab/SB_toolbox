function [mc, M_T0, Phi_T0] = mean_control_cost(A, BtB, mu0, muT, t0, tf, dt)
% M_T0 & Phi_T0 are set to be outputs for further possible computations
n = size(A,1);
stp_n = fix( (tf-t0)/dt ) + 1;

M_t0 = zeros(n,n,stp_n);
M_t_t0 = zeros(n,n,stp_n-1);
for ii = 2:stp_n
    M_t_t0(:,:,ii-1) = A * M_t0(:,:,ii-1) + M_t0(:,:,ii-1)*A' + BtB;
    M_t0(:,:,ii) = M_t0(:,:,ii-1) + M_t_t0(:,:,ii-1)*dt;

end
M_T0 = M_t0(:,:,end);    
    
Phi_Ts_mat = Phi_Ts(A, t0, tf, dt);
Phi_t0_mat = zeros(n,n,stp_n);
for kk = 1:stp_n
    Phi_t0_mat(:,:,kk) = inv(Phi_Ts_mat(:,:,kk))*Phi_Ts_mat(:,:,1);
end
Phi_T0 = Phi_t0_mat(:,:,end);

m_tmp = (muT - Phi_T0 * mu0);
mc = m_tmp' * inv(M_T0) * m_tmp;
end