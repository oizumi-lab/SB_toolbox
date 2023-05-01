function [vtv_t_diag, mean_t_cont_roi, cov_t_cont_roi] = vtv_t_diag_elems(A, BtB, t0, tf, delta_t, Pi_t, m_t, mu_t, cov_t)
% take only diagonal elements of u(t) * u(t)'
% output: 2d array of time steps * dim
    n = size(A,1);
    stp_n = fix( (tf-t0)/delta_t ) + 1;
    
    vtv_t_diag = zeros(stp_n, n);
    mean_t_cont_roi = zeros(stp_n, n);
    cov_t_cont_roi = zeros(stp_n, n);
    for ii = 1:stp_n
        Pi  = Pi_t(:,:,ii);
        cov = cov_t(:,:,ii);
        mu  = mu_t(ii,:)';
        m   = m_t(ii,:)';
        
        term1 = diag( BtB'*Pi*cov*Pi*BtB );
        term2 = diag( BtB'*Pi*(mu*mu')*Pi*BtB );
        term3 = diag( BtB'*Pi*mu*m'*BtB );
        term4 = diag( BtB'*m*m'*BtB );
        
        vtv_t_diag(ii,:) = ( term1 + term2 - 2*term3 + term4)';
        mean_t_cont_roi(ii,:) = (term2-2*term3 + term4)';
        cov_t_cont_roi(ii,:) = term1';
    end


end