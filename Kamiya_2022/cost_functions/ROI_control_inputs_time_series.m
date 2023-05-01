
k = 8; % Task 'Social'
%itr = 1;

mean_ROI_ts_btstrp  = zeros(floor(T/dt)+1, dim, N_btstrp);
cov_ROI_ts_btstrp   = zeros(floor(T/dt)+1, dim, N_btstrp);
total_ROI_ts_btstrp = zeros(floor(T/dt)+1, dim, N_btstrp);

parfor itr = 1:N_btstrp
    fprintf('Bootstrap %d \n', itr)
    
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
    
    mean_ROI_ts_btstrp(:,:,itr) = mean_t_cont_roi_vv;
    cov_ROI_ts_btstrp(:,:,itr)  = cov_t_cont_roi_vv;
    total_ROI_ts_btstrp(:,:,itr)= vtv_t_diag;

end

mean_ROI_ts = mean(mean_ROI_ts_btstrp,3);
cov_ROI_ts  = mean(cov_ROI_ts_btstrp,3);
total_ROI_ts = mean(total_ROI_ts_btstrp,3);

%% Figures
figure;
sgtitle(['Task:', task_list{k}])
subplot(3,1,1)
imagesc(mean_ROI_ts')
title('mean')
colorbar

subplot(3,1,2)
imagesc(cov_ROI_ts')
title('covariance')
colorbar

subplot(3,1,3)
imagesc(total_ROI_ts')
title('total')
colorbar


%% Save figures
times_snapshots = [0:0.1:1];
steps_snapshots = floor(times_snapshots/dt)+1;

niifile_save_dir_name = 'ROI_control_inputs_time_series';
mkdir(niifile_save_dir_name)
for jj = 1:length(steps_snapshots)
    save_nii_file_fixed([niifile_save_dir_name, '/',task_list{k},'_total_step_', num2str(steps_snapshots(jj))], total_ROI_ts(steps_snapshots(jj),:), dim);
    save_nii_file_fixed([niifile_save_dir_name, '/',task_list{k},'_mean_step_',  num2str(steps_snapshots(jj))], mean_ROI_ts(steps_snapshots(jj),:), dim);
    save_nii_file_fixed([niifile_save_dir_name, '/',task_list{k},'_cov_step_',   num2str(steps_snapshots(jj))], cov_ROI_ts(steps_snapshots(jj),:), dim);
    
end



