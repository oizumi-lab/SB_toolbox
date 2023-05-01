function [int_vtv_matrix, int_vtv_mean_matrix, int_vtv_cov_matrix] = ROI_control_distribution(params, emp, bts, P)
%% Control of each ROI (vv)

T  = P.T;
dt = P.dt;

dim      = P.dim;
nTasks   = P.nTasks;
N_btstrp = bts.N_btstrp;

int_vtv_matrix      = zeros(dim, N_btstrp, nTasks);
int_vtv_mean_matrix = zeros(dim, N_btstrp, nTasks);
int_vtv_cov_matrix  = zeros(dim, N_btstrp, nTasks);

parfor itr = 1:N_btstrp
    fprintf('Bootstrap %d \n', itr)
    for k = 1:nTasks
        fprintf('    Task: No. %d \n', k)
        cov1      = emp.cov_rest1(:,:,itr);
        A_est     = params.A_est_matrix(:,:,itr);
        Sigma_est = params.Sigma_est_matrix(:,:,itr);
        cov_task  = emp.cov_tasks(:,:,k,itr);
        mu_task   = emp.mu_tasks(:,k,itr);

        Pi0     = Pi0m_BtB(A_est, Sigma_est, cov1, cov_task, 0, T, dt);
        Pi_t    = Pit_tensor_Pi0_BtB(A_est, Sigma_est, 0, T, dt, Pi0);
        m_t_vec = m_t_BtB(A_est, Sigma_est, zeros(dim,1), mu_task, 0, T, dt, Pi_t);
        mu_t    = mu_t_BtB(A_est, Sigma_est, 0, T, dt, Pi_t, m_t_vec);
        cov_t   = cov_t_BtB(A_est, Sigma_est, cov1, 0, T, dt, Pi_t);

        
        [vtv_t_diag, mean_t_cont_roi_vv, cov_t_cont_roi_vv] =...
            vtv_t_diag_elems(A_est, Sigma_est, 0, T, dt, Pi_t, m_t_vec, mu_t, cov_t);
        
        int_vtv_diag = ( trapz(vtv_t_diag) * dt )';
        int_vtv_matrix(:,itr,k) = int_vtv_diag;
        
        int_vtv_mean_diag = ( trapz(mean_t_cont_roi_vv) * dt )';
        int_vtv_mean_matrix(:,itr,k) = int_vtv_mean_diag;
        
        int_vtv_cov_diag = ( trapz(cov_t_cont_roi_vv) * dt )';
        int_vtv_cov_matrix(:,itr,k) = int_vtv_cov_diag;
        
        
    end
end

%%
%save('~/SBP_matlab/HCP/HCP_bootstrap_concatenate_mat/int_vtv_indiv.mat','int_vtv_matrix', 'int_vtv_mean_matrix','int_vtv_cov_matrix')
end
