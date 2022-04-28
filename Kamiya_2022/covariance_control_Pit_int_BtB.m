function ret = covariance_control_Pit_int_BtB(AA, BtB, Sigma_0, Sigma_T, t0, tf, delta_t, Pit_tensor)
    %%%% int (Q^{-1} BB^{T}) dt
    n_steps = fix( (tf-t0)/delta_t ) + 1;
    d = size(AA,1);
    
    Intgr_vec = zeros(n_steps,1);
    for ii = 1:n_steps
        Intgr_vec(ii,1) = trace( Pit_tensor(:,:,ii)*BtB );
    end
    Intgr = trapz(Intgr_vec,1) * delta_t;
    
    Pi_0 = Pit_tensor(:,:,1);
    Pi_T = Pit_tensor(:,:,end);
    ret = Intgr - trace(Pi_T*Sigma_T-Pi_0*Sigma_0);
end

