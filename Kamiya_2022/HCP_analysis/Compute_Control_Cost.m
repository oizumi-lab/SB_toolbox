function mvt_btstrp = Compute_Control_Cost(P, emp, params)
% P: basic parameters
% -- T
% -- dt
% -- TR
% -- dim
% -- nTasks
% -- N_subj
% emp: structure that contains empirical means & covariances of rest & all
% task states:
% -- mu_tasks
% -- cov_tasks
% params: structure that contains estimated Ornstein Uhlenbeck parameters:
% -- A_est_matrix
% -- Sigma_est_matrix

nTasks   = P.nTasks;                        % Number of tasks
N_btstrp = size(params.A_est_matrix,3);     % Number of bootstrap
dim      = P.dim;                           % Dimension of time series = Number of brain regions

m_btstrp = zeros(nTasks, N_btstrp);         % Store mean control cost values
c_btstrp = zeros(nTasks, N_btstrp);         % Store covariance control cost values

parfor k = 1:nTasks
    for itr = 1:N_btstrp
        %% Assign parameters to compute control costs
        A      = params.A_est_matrix(:,:,itr);       % Drift term of Ornstein Uhlenbeck process
        BtB    = params.Sigma_est_matrix(:,:,itr);   % Diffusion term of OU process
        mu0    = zeros(dim,1);                       % Mean of rest state
        muT    = emp.mu_tasks(:,k,itr);              % Mean of task state
        Sigma0 = emp.cov_tasks(:,:,1,itr);           % Covariance of resting state
        SigmaT = emp.cov_tasks(:,:,k,itr);           % Covariance of task state
        
        %% Compute mean control cost
        fprintf('Computing mean control cost of task no. %d, in bootstrap: %d / %d \n',k, itr, N_btstrp)
        [mc,~,~] = mean_control_cost(A, BtB, mu0, muT, 0, P.T, P.dt);
        m_btstrp(k,itr)   = mc;
        
        %% Compute covariance control cost
        fprintf('Computing covariance control cost of task no. %d, in bootstrap: %d / %d \n',k, itr, N_btstrp)
        % To compute covaraince control cost, we need to perform an
        % integration with a numerically solved matrix dif equation

        Pi0        = Pi0_init_value(A, BtB, Sigma0, SigmaT, 0, P.T, P.dt);    % Initial value (Pi_0) of the differential eq
        Pit_tensor = Pit_numerical(A, BtB, 0, P.T, P.dt, Pi0);                % Numerically solved series of matrices (Pi_t) differential eq, helps computation of covariance control cost
        vc         = covariance_control_cost(BtB, Sigma0, SigmaT, 0, P.T, P.dt, Pit_tensor);
        c_btstrp(k, itr)   = vc;


    end
end

m_btstrp_cct = reshape(m_btstrp, [size(m_btstrp,1),1,size(m_btstrp,2)]);
c_btstrp_cct = reshape(c_btstrp, [size(c_btstrp,1),1,size(c_btstrp,2)]);
mvt_btstrp = cat(2, m_btstrp_cct, c_btstrp_cct, m_btstrp_cct + c_btstrp_cct);


end
