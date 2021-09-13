function [mincost,P] = solveSBP(ini_dist,tar_dist,Q)
%
% Function for computing the minimum control cost and the optimally controlled
% joint end point distribution P
%
% Input: 
%   - ini_dist : initial distribution (k x 1 vector)
%   - tar_dist : target distribution (k x 1 vector) 
%   - Q : uncontrolled endpoint distribution (k x k matrix)
% 
% Output: 
%   - KLD: minimum control cost, KL divergence between the uncontrolled 
%          and optimally controlled distributions.
%   - P : optimally controlled endpoint distribution (k x k matrix)
%   
% Note that k is the number of states
%

%% Finding the optimally controlled endpoint distribution
k = size(ini_dist,1);
U = Q.*(-log(Q));
lambda = 1;
[D,L,u,v]=sinkhornTransport(ini_dist,tar_dist,Q,U,lambda,'marginalDifference');

P = bsxfun(@times,v',(bsxfun(@times,u,Q)));
P = P/sum(sum(P));

%% Computing the minimum control cost
mincost = compute_KLD(P(:),Q(:));

end