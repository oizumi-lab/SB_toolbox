function save_vtv_ROI_vec_to_nii(niifile_save_dir_name, int_vtv_matrix, int_vtv_mean_matrix, int_vtv_cov_matrix)
%% Save nii file for Connectome Workbench visualitzation
all_task_file_list =...
    {'REST1','REST2','EMOTION', 'GAMBLING','LANGUAGE', 'MOTOR',...
    'RELATIONAL','SOCIAL','WM_0bk','WM_2bk','WM'};
task_file_name_detailed = {'REST1', 'REST2', 'EMOT_fear', 'EMOT_neut',...
    'GAMB_win', 'GAMB_loss', 'LANG_story', 'LANG_math',...
    'MOT_lf', 'MOT_rf', 'MOT_lh', 'MOT_rh', 'MOT_t',...
    'RELAT_relation', 'RELAT_match', 'SOC_mental', 'SOC_rnd',...
    'WM0bk_body', 'WM0bk_faces', 'WM0bk_places', 'WM0bk_tools',...
    'WM2bk_body', 'WM2bk_faces', 'WM2bk_places', 'WM2bk_tools'};

nTasks = size(int_vtv_matrix,3);
dim    = size(int_vtv_matrix,1);
total_vtv_ave = squeeze( mean(int_vtv_matrix,2) );
mean_vtv_ave  = squeeze( mean(int_vtv_mean_matrix,2) );
cov_vtv_ave   = squeeze( mean(int_vtv_cov_matrix,2) );

if nTasks == 25
    task_file_name = task_file_name_detailed;
elseif nTasks <= 11
    task_file_name = all_task_file_list;
end

for k = 1:nTasks-2
    task_vtv_total_array = total_vtv_ave(:,k+2);
    task_vtv_mean_array  = mean_vtv_ave(:,k+2);
    task_vtv_cov_array   = cov_vtv_ave(:,k+2);

    mkdir(niifile_save_dir_name)
    file_nii_name_t = strcat(niifile_save_dir_name, '/HCP_ROImap_detailed_total_v_',task_file_name(k+2));
    file_nii_name_t = file_nii_name_t{1,1};
    file_nii_name_m = strcat(niifile_save_dir_name, '/HCP_ROImap_detailed_mean_v_',task_file_name(k+2));
    file_nii_name_m = file_nii_name_m{1,1};
    file_nii_name_c = strcat(niifile_save_dir_name, '/HCP_ROImap_detailed_cov_v_',task_file_name(k+2));
    file_nii_name_c = file_nii_name_c{1,1};
    
    
    save_nii_file_fixed(file_nii_name_t, task_vtv_total_array, dim);
    save_nii_file_fixed(file_nii_name_m, task_vtv_mean_array, dim);
    save_nii_file_fixed(file_nii_name_c, task_vtv_cov_array, dim);
    
end

end

