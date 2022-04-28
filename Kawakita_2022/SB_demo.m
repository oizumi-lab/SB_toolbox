%
% Demo program for computing brain state transition cost from
% the resting state to 0-back and 2-back tasks
% This code generates Fig. 3 (a) in
% Quantifying brain state transition cost via Schrödinger’s bridge
% Genji Kawakita, Shunsuke Kamiya, Shuntaro Sasai, Jun Kitazono, Masafumi Oizumi
% bioRxiv: https://doi.org/10.1101/2021.05.24.445394
% 

%% 1. Loading data
load('./demo_data.mat');
Q = Data.Q; % the uncontrolled path (the joint dsitribution, q(X_0,X_T))
rest_dist = Data.probDists(1,:); % probability distribution for the rest
wm0_dist = Data.probDists(8,:); % probabilty distribution for 0-back task
wm2_dist = Data.probDists(9,:); % probability distribution for 2-back task

%% 2. Computing brain state transition cost from the resting state to 
%     0-back and 2-back tasks

% cost for rest -> 0-back
cost_0back = solveSBP(rest_dist',wm0_dist',Q);
% cost for rest -> 2-back
cost_2back = solveSBP(rest_dist',wm2_dist',Q);

%% 3. Plotting the result

col_blue = [114 147 203]/255;
col_red = [211 94 96]/255;

figure
z = categorical({'Rest -> 0-back' 'Rest -> 2-back'});
z1 = categorical({'Rest -> 0-back'});
z2 = categorical({'Rest -> 2-back'});
rest2wms = [cost_0back;cost_2back];
b3 = bar(z1,rest2wms(1,1));
set(b3,'FaceColor',col_blue)
hold on
b4 = bar(z2,rest2wms(2,1));
set(b4,'FaceColor',col_red)
hold on

title("Transition cost from the rest to 0-back and 2-back tasks");
ylabel('State Transition Cost','FontSize',16);
axis tight
box off
hold off
