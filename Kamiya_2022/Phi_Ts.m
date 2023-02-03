function Phi_Ts_mat = Phi_Ts(A, t0, tf, dt)
    n = size(A,1);
    stp_n = round( (tf-t0)/dt ) + 1;
    
    Id = eye(n);
    Phi_Ts_mat = zeros(n,n,stp_n);
    Phi_Q_s_Ts = zeros(n,n,stp_n-1);
    Mhat_T0_mat = zeros(n,n,stp_n);
    Phi_Ts_mat(:,:,end) = Id;
    for ii = stp_n:-1:2
        Phi_Q_s_Ts(:,:,ii-1) = - Phi_Ts_mat(:,:,ii) * A;
        Phi_Ts_mat(:,:,ii-1) = Phi_Ts_mat(:,:,ii) - Phi_Q_s_Ts(:,:,ii-1) * dt;
    end
    
end

