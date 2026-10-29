function popchange = opioid(t,pops,p);

S = pops(p.k+1);
A = pops(p.k+2);

popchange = zeros(p.num_classes,1);

popchange(1) = p.rho*S - (p.alpha*p.k + p.sigma_vec(1))*pops(1) - p.beta * pops(1) * A / p.N - p.mu * pops(1);
if p.cutoff == 1
else
    for i = 2:p.cutoff
        popchange(i) = p.alpha*p.k*pops(i-1) - (p.alpha*p.k + p.sigma_vec(i))*pops(i) - p.beta * pops(i) * A / p.N - p.mu * pops(i);
        popchange(p.k + 4 + i) = p.sigma_vec(i)*pops(i) + p.beta * pops(i) * A / p.N;
    end
end
    popchange(p.k + 5) = p.sigma_vec(1)*pops(1) + p.beta * pops(1) * A / p.N;

    popchange(p.k+1) = -p.rho*S + p.alpha*pops(p.cutoff)*p.k - p.beta * S * A / p.N + p.eta*A - p.mu*S + p.mu * (S + sum(pops(1:p.k)) + A) + p.muOD * A;
    popchange(p.k+2) =  sum(pops(1:p.cutoff)'.*p.sigma_vec(1:p.cutoff)) + p.beta * (S + sum(pops(1:p.cutoff))) * A / p.N - p.eta*A - A * (p.mu + p.muOD);

    popchange(p.k+3) = p.rho*S;
    popchange(p.k+4) = sum(pops(1:p.k)'.*p.sigma_vec(1:end))+ + p.beta * (S + sum(pops(1:p.cutoff))) * A / p.N;

    popchange(2*p.k + 5) = sum(pops(1:p.k)'.*p.sigma_vec(1:end));
    popchange(2*p.k + 6) = p.beta*S*A/p.N;
    popchange(2*p.k + 7) = p.beta*sum(pops(1:p.cutoff))*A/p.N;
