function Pi_t = Pit_numerical(A, BtB, t0, tf, dt, Pi0)
    n = size(A,1);
    n_steps = round( (tf-t0)/dt ) + 1;
    Pi_t = zeros(n,n,n_steps);
    
    Pi_t(:,:,1) = Pi0; 
    for ii = 2:n_steps
        k1 = - A'* Pi_t(:,:,ii-1) -  Pi_t(:,:,ii-1)*A...
            + Pi_t(:,:,ii-1)*BtB*Pi_t(:,:,ii-1);
        k2 = - A'*(Pi_t(:,:,ii-1) + dt*k1/2) - (Pi_t(:,:,ii-1) + dt*k1/2)*A...
            + (Pi_t(:,:,ii-1) + dt*k1/2)*BtB*(Pi_t(:,:,ii-1) + dt*k1/2);
        k3 = - A'*(Pi_t(:,:,ii-1) + dt*k2/2) - (Pi_t(:,:,ii-1) + dt*k2/2)*A...
            + (Pi_t(:,:,ii-1) + dt*k2/2)*BtB*(Pi_t(:,:,ii-1) + dt*k2/2);
        k4 = - A'*(Pi_t(:,:,ii-1) + dt*k3)   - (Pi_t(:,:,ii-1) + dt*k3)*A...
            + (Pi_t(:,:,ii-1) + dt*k3)*BtB*(Pi_t(:,:,ii-1) + dt*k3);
        Pi_t(:,:,ii) = Pi_t(:,:,ii-1) + dt/6 * (k1+2*k2+2*k3+k4);
    end
end

