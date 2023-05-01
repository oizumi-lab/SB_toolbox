function Pi_t = Pit_tensor_Pi0_BtB(AA, BtB, t0, tf, delta_t, Pi0)
    n = size(AA,1);
    n_steps = round( (tf-t0)/delta_t ) + 1;
    Pi_t = zeros(n,n,n_steps);
    
    Pi_t(:,:,1) = Pi0; 
    for ii = 2:n_steps
        k1 = - AA'* Pi_t(:,:,ii-1) -  Pi_t(:,:,ii-1)*AA...
            + Pi_t(:,:,ii-1)*BtB*Pi_t(:,:,ii-1);
        k2 = - AA'*(Pi_t(:,:,ii-1) + delta_t*k1/2) - (Pi_t(:,:,ii-1) + delta_t*k1/2)*AA...
            + (Pi_t(:,:,ii-1) + delta_t*k1/2)*BtB*(Pi_t(:,:,ii-1) + delta_t*k1/2);
        k3 = - AA'*(Pi_t(:,:,ii-1) + delta_t*k2/2) - (Pi_t(:,:,ii-1) + delta_t*k2/2)*AA...
            + (Pi_t(:,:,ii-1) + delta_t*k2/2)*BtB*(Pi_t(:,:,ii-1) + delta_t*k2/2);
        k4 = - AA'*(Pi_t(:,:,ii-1) + delta_t*k3)   - (Pi_t(:,:,ii-1) + delta_t*k3)*AA...
            + (Pi_t(:,:,ii-1) + delta_t*k3)*BtB*(Pi_t(:,:,ii-1) + delta_t*k3);
        Pi_t(:,:,ii) = Pi_t(:,:,ii-1) + delta_t/6 * (k1+2*k2+2*k3+k4);
    end
end

