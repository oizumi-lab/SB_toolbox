function vc = covariance_control_cost(BtB, Sigma0, SigmaT, t0, tf, dt, Pit_tensor)
    %%%% int (Q^{-1} BB^{T}) dt
    n_steps = fix( (tf-t0)/dt ) + 1;
    
    Intgr_vec = zeros(n_steps,1);
    for ii = 1:n_steps
        Intgr_vec(ii,1) = trace( Pit_tensor(:,:,ii)*BtB );
    end
    Intgr = trapz(Intgr_vec,1) * dt;
    
    Pi_0 = Pit_tensor(:,:,1);
    Pi_T = Pit_tensor(:,:,end);
    vc = Intgr - trace(Pi_T*SigmaT-Pi_0*Sigma0);
end

