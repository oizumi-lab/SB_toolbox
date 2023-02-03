% Toy example code. Compute mean & covariance control costs from "state 1"
% to "state 2".
%   - State 1: initial brain state, stationary multivariate Ornstein Uhlenbeck process
%   - State 2: target brain state

%% Step 0: simulate time serireses ("state 1" and "state 2")
% "state1": a two-dimensonal stationary time series simulated
% according to a multivariate OU equation:
%           dx = A x dt + B dw,
% where A, B \in R^{2 x 2}. 
% "state2": arbitrary time series

rng(9876)

dim     = 2;                    % state dimension
A       = [1 2; -3 -4];         % Drift term
BtB     = [3 1 ;1 1];           % Defines diffusion matrix B (BtB := B * B')

Sigma_0 = lyap(A, BtB);         % Stationary covariance of state1, see e.g. "Applied Stochastic Differential Equations" by Sarkka & Arno for further derivation
x_0     = 0.1 * randn(1,dim);   % Initial value of state 1

dt      = 1e-3;                 % Time interval to simulate MOU process & compute integration
TR      = 1e-2;                 % Time interval to sample time series
T       = 100;                  % Time length of time series
N       = floor(T/dt) + 1;      % # of time steps determined by T and dt

% simulate state 1 (stationary OU process)
state1       = simulate_MOU(x_0, A, BtB, TR, dt, T);

% simulate state 2 
mu_state2    = [2 -1];
Sigma_state2 = [1 0;0 2];
state2       = repmat(mu_state2,N,1) + randn(N,dim) * (Sigma_state2)^(1/2);

%% Step 1: compute empirical mean & covariance of each state
% As state1.mat is sumulated by a stationary time series, the mean is
% strictly zero. Thus the empirical mean also should be close to zero.

mu_1    = mean(state1,1)';
Sigma_1 = cov(state1);
mu_2    = mean(state2,1)';
Sigma_2 = cov(state2);

%% Step 2: estimate Ornstein Uhlenbeck parameters
% You only need OU parameters of the uncontrolled state (here, state 1) to
% compute control cost. Here we employ lasso regression.
lambda =  0.0001;           % lasso hyperparameter. One can use CV for deciding this value
Y = state1(2:end,:);
X = state1(1:end-1,:);

A_d_hat = zeros(dim,dim);   % This should estimate "expm(A*TR)"
for kk = 1:dim
    A_kk          = lasso(X,Y(:,kk),'Lambda',lambda);
    A_d_hat(kk,:) = A_kk';
end

A_hat = logm(A_d_hat)/TR;   % Estimation of A
BtB_hat = -( A_hat * Sigma_1 + Sigma_1 * A_hat'); % Estimation of BtB

%% Step 3: compute mean & covariance control costs
[mean_cost,~,~] = mean_control_cost(A_hat, BtB_hat, mu_1, mu_2, 0, 1, dt);

Pi0        = Pi0_init_value(A_hat, BtB_hat, Sigma_1, Sigma_2, 0, 1, dt);    % Initial value (Pi_0) of the differential eq
Pit_tensor = Pit_numerical(A, BtB, 0, 1, dt, Pi0);  
cov_cost   = covariance_control_cost(BtB, Sigma_1, Sigma_2, 0, 1, dt, Pit_tensor);
        
total_cost = mean_cost + cov_cost;

fprintf(' ---------------------------------------- \n')
fprintf(' Mean control cost       = %5.3d, \n', mean_cost)
fprintf(' Covariance control cost = %5.3d, \n', cov_cost)
fprintf(' Total control cost      = %5.3d. \n', total_cost)

