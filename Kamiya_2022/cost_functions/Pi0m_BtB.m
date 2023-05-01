function ret = Pi0m_BtB(AA, BtB, Sigma_0, Sigma_T, t0, tf, delta_t)
    n_steps = fix( (tf-t0)/delta_t ) + 1;
    d = size(AA,1);

    N_T0_mat = zeros(d,d,n_steps);
    for ii = 1:n_steps
        s_i = t0 + (ii-1) * delta_t;
        EXP = expm(AA*(t0-s_i));
        N_T0_mat(:,:,ii) = EXP * BtB * EXP';
    end
    N_T0 = trapz(N_T0_mat, 3) * delta_t;
    
    Ninvsq =  N_T0^(-1/2); EXP_tft0 = expm(-AA*(tf-t0));
    S0 = Ninvsq * Sigma_0 * Ninvsq;
    ST = Ninvsq * EXP_tft0 * Sigma_T * EXP_tft0' * Ninvsq;
    Sinvsq = S0^(-1/2);
    %INV_SN = Ninvsqrt * S0^(-1/2);
    
    %ret = INV_SN * (S0 + (1/2)*eye(d) - (S0^(1/2)*ST*S0^(1/2) + (1/4)*eye(d))^(1/2)) * INV_SN;
    ret =  Ninvsq * Sinvsq * (S0 + (1/2)*eye(d) - (S0^(1/2)*ST*S0^(1/2) + (1/4)*eye(d))^(1/2)) * Sinvsq * Ninvsq;
end
