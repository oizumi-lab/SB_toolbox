function mvt_btstrp = Compute_Control_Cost(P, emp, params)
% emp: contains
% - mu_tasks
% - cov_tasks
% params: contains
% - A_est_matrix
% - Sigma_est_matrix

nTasks   = size(emp.mu_tasks,2);
N_btstrp = size(params.A_est_matrix,3);
dim      = size(params.A_est_matrix,1);

m_btstrp = zeros(nTasks, N_btstrp);
c_btstrp = zeros(nTasks, N_btstrp);

parfor k = 1:nTasks
    for itr = 1:N_btstrp
        fprintf('Computing mean control cost of task no. %d, in bootstrap: %d / %d \n',k, itr, N_btstrp)
        A      = params.A_est_matrix(:,:,itr);
        BtB    = params.Sigma_est_matrix(:,:,itr);
        Sigma0 = emp.cov_tasks(:,:,1,itr);
        %Sigma0 = emp.cov_rest1(:,:,itr);
        SigmaT = emp.cov_tasks(:,:,k,itr);
        
        %%
        [mc,~,~] = mean_control_deterministic(A, BtB, zeros(dim,1), emp.mu_tasks(:,k,itr), 0, P.T, P.dt);
        m_btstrp(k,itr)   = mc;
        
        %%
        fprintf('Computing covariance control cost of task no. %d, in bootstrap: %d / %d \n',k, itr, N_btstrp)
        Pi0 = Pi0m_BtB(A, BtB, Sigma0, SigmaT, 0, P.T, P.dt);
        Pit_tensor = Pit_tensor_Pi0_BtB(A, BtB, 0, P.T, P.dt, Pi0);
        vc = covariance_control_Pit_int_BtB(A, BtB, Sigma0, SigmaT, 0, P.T, P.dt, Pit_tensor);
        c_btstrp(k, itr)   = vc;


    end
end

m_btstrp_cct = reshape(m_btstrp, [size(m_btstrp,1),1,size(m_btstrp,2)]);
c_btstrp_cct = reshape(c_btstrp, [size(c_btstrp,1),1,size(c_btstrp,2)]);
mvt_btstrp = cat(2, m_btstrp_cct, c_btstrp_cct, m_btstrp_cct + c_btstrp_cct);


end
