function ret = Pi0_init_value(A, BtB, Sigma0, SigmaT, t0, tf, dt)
    n_steps = fix( (tf-t0)/dt ) + 1;
    d = size(A,1);

    N_T0_mat = zeros(d,d,n_steps);
    for ii = 1:n_steps
        s_i = t0 + (ii-1) * dt;
        EXP = expm(A*(t0-s_i));
        N_T0_mat(:,:,ii) = EXP * BtB * EXP';
    end
    N_T0 = trapz(N_T0_mat, 3) * dt;
    
    Ninvsq =  N_T0^(-1/2); EXP_tft0 = expm(-A*(tf-t0));
    S0 = Ninvsq * Sigma0 * Ninvsq;
    ST = Ninvsq * EXP_tft0 * SigmaT * EXP_tft0' * Ninvsq;
    Sinvsq = S0^(-1/2);
    %INV_SN = Ninvsqrt * S0^(-1/2);
    
    %ret = INV_SN * (S0 + (1/2)*eye(d) - (S0^(1/2)*ST*S0^(1/2) + (1/4)*eye(d))^(1/2)) * INV_SN;
    ret =  Ninvsq * Sinvsq * (S0 + (1/2)*eye(d) - (S0^(1/2)*ST*S0^(1/2) + (1/4)*eye(d))^(1/2)) * Sinvsq * Ninvsq;
end
