function change = opioid_basic(~,states,p)
P = states(1:p.k);
S = states(p.k+1);
A = states(p.k+2);
social = p.beta_const + p.beta*A/p.N;
hazard = p.sigma_vec(:) + social;
change = zeros(p.num_classes,1);

change(1) = p.rho*S - (p.alpha*p.k + hazard(1) + p.mu)*P(1);
for i = 2:p.cutoff
    change(i) = p.alpha*p.k*P(i-1) - ...
        (p.alpha*p.k + hazard(i) + p.mu)*P(i);
end

new_addictions = sum(P(1:p.cutoff).*hazard(1:p.cutoff)) + social*S;
population = S + sum(P) + A;
change(p.k+1) = -(p.rho+social)*S + p.alpha*p.k*P(p.cutoff) + ...
    p.eta*A - p.mu*S + p.mu*population + p.muOD*A;
change(p.k+2) = new_addictions - p.eta*A - (p.mu+p.muOD)*A;
change(p.k+3) = p.rho*S;
change(p.k+4) = new_addictions;
end
