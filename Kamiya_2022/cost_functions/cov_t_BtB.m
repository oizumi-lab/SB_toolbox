function cov_t = cov_t_BtB(A, BtB, cov_0, t0, tf, delta_t, Pi_t)

    n = size(A,1);
    stp_n = fix( (tf-t0)/delta_t ) + 1;
    
    cov_t = zeros(n, n, stp_n);
    cov_t(:,:,1) = cov_0;
    for ii = 2:stp_n
        F = A - BtB * Pi_t(:,:,ii-1);
        k1 = F *  cov_t(:,:,ii-1) +  cov_t(:,:,ii-1) * F' + BtB;
        k2 = F * (cov_t(:,:,ii-1) + delta_t*k1/2) + (cov_t(:,:,ii-1) + delta_t*k1/2) * F' + BtB;
        k3 = F * (cov_t(:,:,ii-1) + delta_t*k2/2) + (cov_t(:,:,ii-1) + delta_t*k2/2) * F' + BtB;
        k4 = F * (cov_t(:,:,ii-1) + delta_t*k3)   + (cov_t(:,:,ii-1) + delta_t*k3) * F' + BtB;
        cov_t(:,:,ii) = cov_t(:,:,ii-1) + delta_t/6 * (k1+2*k2+2*k3+k4);
    end
    
end