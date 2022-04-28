function KLD = compute_KLD(P,Q)

% Function for computing the Kullback-Leibler divergence between P and Q

KLD = 0;
for i=1:length(P)
    if P(i) ~= 0
        KLD = KLD + P(i)*(log(P(i))-log(Q(i)));
    end
end

end
