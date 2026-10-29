function popchange = opioid_relapse(t,pops,p)

S = pops(p.k+1);
A = pops(p.k+2);
R = pops(p.k+3);

popchange = zeros(p.num_classes,1);

%% Prescription classes
popchange(1) = p.rho*S - (p.alpha*p.k + p.sigma_vec(1)+p.beta_const+p.beta*A/p.N)*pops(1) - p.mu * pops(1);

if p.cutoff > 1
    for i = 2:p.cutoff
        popchange(i) = p.alpha*p.k*pops(i-1) ...
                     - (p.alpha*p.k + p.sigma_vec(i)+p.beta_const+p.beta*A/p.N)*pops(i) - p.mu * pops(i);
        popchange(p.k + 4 + i) = (p.sigma_vec(i)+p.beta_const+p.beta*A/p.N)*pops(i);
    end
end

popchange(p.k + 5) = (p.sigma_vec(1)+p.beta_const+p.beta*A/p.N)*pops(1);

%% Compartments

popchange(p.k+1) = -(p.rho + p.beta_const+p.beta*A/p.N)*S + p.alpha*pops(p.cutoff)*p.k - p.mu*S + p.mu * (S + sum(pops(1:p.k)) + A + R) + p.muOD * A;

popchange(p.k+2) = sum(pops(1:p.cutoff)'.*(p.sigma_vec(1:p.cutoff)+p.beta_const+p.beta*A/p.N)) ...
                 + (p.beta_const+p.beta*A/p.N)*S ...
                 - p.eta*A ...
                 + p.delta*R ...
                 - A * (p.mu + p.muOD);

popchange(p.k+3) = p.eta*A - p.delta*R - p.mu*R;

popchange(2*p.k+4) = p.rho*S;

popchange(2*p.k+5) = sum(pops(1:p.k)'.*(p.sigma_vec(1:end)+p.beta_const+p.beta*A/p.N)) + (p.beta_const+p.beta*A/p.N)*S;

%% Cumulative relapse addictions
popchange(2*p.k+6) = p.delta*R;

popchange(2*p.k + 7) = sum(pops(1:p.k)'.*p.sigma_vec(1:end));

popchange(2*p.k + 8) = (p.beta_const+p.beta*A/p.N)*S + (p.beta_const+p.beta*A/p.N)*sum(pops(1:p.k));

popchange(2*p.k + 9) = (p.beta_const + p.beta*A/p.N) * S;

popchange(2*p.k + 10) = (p.beta_const + p.beta*A/p.N) * sum(pops(1:p.k));

end
