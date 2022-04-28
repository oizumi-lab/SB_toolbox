%% Basic settings

nTasks = 11;
%HCP_load_vars_new

T   = 1.0;
dt  = 1e-3;
TR  = 0.72;
dim = 100;
N_subj = 352;

task_list = {'REST1', 'REST2', 'EMOTION', 'GAMBLING', 'LANGUAGE',...
                'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM'};

%% compute mean and covariance
P.T   = T;
P.dt  = dt;
P.TR  = TR;
P.dim = dim;
P.nTasks = nTasks;
%P.database_name_cell = database_name_cell;
P.N_subj = N_subj;

%% Bootstrapping samples
N_samples = 100;
N_btstrp  = 100;

rng(1234);
M = randi(N_subj, N_samples, N_btstrp); % sample arrays are column vectors

bts.N_samples = N_samples;
bts.N_btstrp  = N_btstrp;
bts.M         = M;


%% load data
load('Kamiya_2022/Data.mat')


%% compute control energy
mvt_btstrp = Compute_Control_Cost(P, emp, params);

%% show bar figures
FigPars.tasks_included_vec = [3:8,11];          % change this to select tasks to include in graph
FigPars.task_detailed = double( nTasks >= 25 ); % choose 0 or 1 to change categories of tasks
FigPars.bars_to_show = 'only_total';            % choose from 'all', 'only_total', 'only_mean', 'only_covariance'
FigPars.bars_order   = 'ascending';             % choose from 'normal' or 'ascending'
close all
[f,lgd] = HCP_bootstrap_concatenate_makebargraph(FigPars, nTasks, mvt_btstrp);
title([FigPars.bars_to_show(6:end), ' control costs'])

set(gca, 'FontSize',18);
set(lgd, 'FontSize',18);

%% Compute distributions
[int_vtv_matrix, int_vtv_mean_matrix, int_vtv_cov_matrix] = ROI_control_distribution(params, emp, bts, P);
% int_vtv_matrix: dim * N_btstrp * nTasks 

int_vtv_matrix      = 1e5 * int_vtv_matrix;
int_vtv_mean_matrix = 1e5 * int_vtv_mean_matrix;
int_vtv_cov_matrix  = 1e5 * int_vtv_cov_matrix;

niifile_save_dir_name = 'Kamiya_2022/HCP_ROI_control_map';
save_vtv_ROI_vec_to_nii(niifile_save_dir_name, int_vtv_matrix, int_vtv_mean_matrix, int_vtv_cov_matrix)

%% Preparation before computation

ratio_selection = 0.3;
rank_ind = 10;
tasks_to_choose = [3:8,11];
vtv_cell       = {int_vtv_mean_matrix(:,:,tasks_to_choose), int_vtv_cov_matrix(:,:,tasks_to_choose), int_vtv_matrix(:,:,tasks_to_choose)};
file_name_cell = {'mean', 'covariance', 'total'};

%% Common ROIs important for input
for ll = 1:3
    vtv_matrix         = vtv_cell{ll};
    file_name_percents = [niifile_save_dir_name, '/ROImap_', file_name_cell{ll} ,'_top30percent_times_appear_seven_tasks'];
    file_name_topROIs  = [niifile_save_dir_name, '/ROImap_', file_name_cell{ll} ,'_top_rank_10_seven_tasks'];
    file_name_meantasks= [niifile_save_dir_name, '/ROImap_', file_name_cell{ll} ,'_simple_task_means_seven_tasks'];
    
    ROIs_top_percents_to_control(vtv_matrix, ratio_selection, file_name_percents);
    ROIs_top_howmany_to_control( vtv_matrix, rank_ind,        file_name_topROIs)
    ROIs_tasks_mean_control(     vtv_matrix, file_name_meantasks)
end

%% top ten ROI's - put colors on ROI's only over four
for ll = 1:3
    vtv_matrix         = vtv_cell{ll};
    file_name_topROIs_over4  = [niifile_save_dir_name, '/ROImap_', file_name_cell{ll} ,'_top_rank_10_seven_tasks_over4'];
    ROIs_top_howmany_to_control( vtv_matrix, rank_ind, file_name_topROIs_over4, 4)
end


%% Select significant ROIs in all tasks average
mean_int_vtv_all_tasks = cell(3,1);
average_rank_ROIs = cell(3,1);
randk_ind = 10;

for ll = 1:3
    vtv_matrix         = vtv_cell{ll};
    
    %% select fixed number of top ROIs
    average_vec = mean( squeeze( mean(vtv_matrix,2) ),2);
    average_numbered_vec = [[1:dim]', mean( squeeze( mean(vtv_matrix,2) ),2)];
    mean_int_vtv_all_tasks{ll} = average_numbered_vec;
   
    average_rank_vec = sortrows(average_numbered_vec, 2, 'descend');
    average_rank_ROIs{ll} = average_rank_vec(1:rank_ind,1);
    
    %% select top ROIs that explains the fixed number of percent of the whole input
    ROI_control = [ [1:dim]', average_vec];
    sum_all_ROI = sum(ROI_control(:,2));
    ratio_sum   = ratio_selection * sum_all_ROI;
    
    ROI_control_sorted     = sortrows(ROI_control, 2, 'descend');
    cumulative_sorted_vec = cumsum(ROI_control_sorted(:,2));
    ind = min( find( cumulative_sorted_vec >= ratio_sum ));
    selected_ROIs = ROI_control_sorted(1:ind,1);
    
    select_vec = zeros(dim,1);
    select_vec(selected_ROIs) = 1;
    
    average_vec_only_selected_ROIs = average_vec .* select_vec;
    
    save_nii_file_fixed([niifile_save_dir_name, '/ROImap_', file_name_cell{ll} ,'_top',num2str(ratio_selection*100),'percent_of_average_seven_tasks'], average_vec_only_selected_ROIs, dim);
    
    %main_control_ROIs_cell{k,1} = selected_ROIs;
end