%% Basic settings

addpath(genpath("../"))

T      = 1.0;      % Final time (sec)
dt     = 1e-3;     % Time interval to compute integrations
TR     = 0.72;     % TR of fMRI data
dim    = 100;      % Number of data dimensions( = # of brain regions)
N_subj = 352;      % Number of subjects contained

task_list = {'REST1', 'REST2', 'EMOTION', 'GAMBLING', 'LANGUAGE',...
                'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM'};
nTasks = length(task_list); % Number of tasks

%% Assign basic parameters into structure P
P.T   = T;
P.dt  = dt;
P.TR  = TR;
P.dim = dim;
P.nTasks = nTasks;
P.N_subj = N_subj;

%% Bootstrapping samples
% In this section, we randomly choose "N_samples" subjects for "N_btstrp"
% times to create randomly assigend bootstrap samples. 
% To run this code, we can omit this section as we are given "Data.mat" 
% that contains means & covariances of % the whole randomly assigned 
% subjects.
N_samples = 100;
N_btstrp  = 100;

rng(1234);
M = randi(N_subj, N_samples, N_btstrp); % sample arrays are column vectors

bts.N_samples = N_samples;
bts.N_btstrp  = N_btstrp;
bts.M         = M;


%% load data
load('../Data.mat')

%% compute control energy
mvt_btstrp = Compute_Control_Cost(P, emp, params);

%% show bar figures
tasks_included_vec = [1:nTasks];          % change this to select tasks to include in graph

close all
[f,lgd] = control_cost_make_bargraph(tasks_included_vec, mvt_btstrp);
title('Control costs')

set(gca, 'FontSize',18);
set(lgd, 'FontSize',18);

