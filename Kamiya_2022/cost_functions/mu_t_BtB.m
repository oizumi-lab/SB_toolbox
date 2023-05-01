function mu_t = mu_t_BtB(A, BtB, t0, tf, delta_t, Pi_t, m_t)

    n = size(A,1);
    stp_n = fix( (tf-t0)/delta_t ) + 1;
    
    mu_t = zeros(stp_n,n);
    % mu_t(1,:)' = [0;0;...;0];
    
    for ii = 2:stp_n
        F = A - BtB * Pi_t(:,:,ii-1);
        k1 = F *  mu_t(ii-1,:)'                 + BtB * m_t(ii-1,:)';
        k2 = F * (mu_t(ii-1,:)' + delta_t*k1/2) + BtB * m_t(ii-1,:)';
        k3 = F * (mu_t(ii-1,:)' + delta_t*k2/2) + BtB * m_t(ii-1,:)';
        k4 = F * (mu_t(ii-1,:)' + delta_t*k3)   + BtB * m_t(ii-1,:)';
        mu_t(ii,:) = ( mu_t(ii-1,:)' + delta_t/6 * (k1+2*k2+2*k3+k4) )';
    end
    
end