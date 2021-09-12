function [P,KLD] = compute_cost(ini_dist,tar_dist,Q)
%
% Function for computing the optimally controlled joint endpoint
% distribution and transition cost.
%
% Input: 
%   - ini_dist : initial distribution (k x 1 vector)
%   - tar_dist : target distribution (k x 1 vector) 
%   - Q : uncontrolled endpoint distribution (k x k matrix)
% Output: 
%   - P : optimally controlled endpoint distribution (k x k matrix)
%   - KLD: transition cost, KL divergence between the uncontrolled 
%          and controlled paths.
% Note that k is the number of brain states (clusters)
%

%% Computing P*
k = size(ini_dist,1);
U = Q.*(-log(Q));
lambda = 1;
[D,L,u,v]=sinkhornTransport(ini_dist,tar_dist,Q,U,lambda,'marginalDifference');

P = bsxfun(@times,v',(bsxfun(@times,u,Q)));
P = P/sum(sum(P));

%% Computing KLD
KLD = 0.0;
for i=1:k
    for j=1:k
        KLD = KLD + P(i,j)*(log2(P(i,j))-log2(Q(i,j)));
    end
end

end