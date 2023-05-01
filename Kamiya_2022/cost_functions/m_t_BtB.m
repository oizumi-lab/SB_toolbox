function m_t_vec = m_t_BtB(A, BtB, mu_0, mu_T, t0, tf, delta_t, Pit_tensor)
% This function computes the supplementary function m(t) in Chen et al
% (2016). m(t) is NOT the mean of x(t) at time t !!!! For finding the
% mean, utilize function named mu_t_BtB instead.
% (The final value of m(t) (=m(T)) is identical to mu_T, which is the mean
% of the target distribution.

    n = size(A,1);
    stp_n = fix( (tf-t0)/delta_t ) + 1;
    
    Id = eye(n);
    Phi_Q_Ts    = zeros(n,n,stp_n);
    Phi_Q_s_Ts  = zeros(n,n,stp_n-1);
    Mhat_T0_mat = zeros(n,n,stp_n);
    Phi_Q_Ts(:,:,end) = Id;
    for ii = stp_n:-1:2
        Phi_Q_s_Ts(:,:,ii-1) = - Phi_Q_Ts(:,:,ii) * (A - BtB*Pit_tensor(:,:,ii));
        Phi_Q_Ts(:,:,ii-1) = Phi_Q_Ts(:,:,ii) - Phi_Q_s_Ts(:,:,ii-1) * delta_t;
    end
    Phi_Q_T0 = Phi_Q_Ts(:,:,1);
    
    for ii = 1:stp_n
        Mhat_T0_mat(:,:,ii) = Phi_Q_Ts(:,:,ii)*BtB*Phi_Q_Ts(:,:,ii)';
    end
    Mhat_T0 = trapz(Mhat_T0_mat,3) * delta_t;
    inv_Mhat_T0 = inv(Mhat_T0);
    
    m_t_vec = zeros(stp_n, n);
    for ii = 1:stp_n
        m_t_vec(ii,:) = ( Phi_Q_Ts(:,:,ii)' * inv_Mhat_T0 * (mu_T - Phi_Q_T0 * mu_0) )';
    end
    
end