function popchange = opioid(t,pops,p);


% pops(1) - pops(p.k) are the P classes
S = pops(p.k+1);
A = pops(p.k+2);


popchange = zeros(p.num_classes,1);


% Model equations

if p.k == 1

    % prescription class
    popchange(1) = p.rho*S - p.alpha*pops(1) - p.sigma_avg*pops(1); 

    % suceptible, addicted
    popchange(2) = - p.rho*S + p.alpha*pops(1);
    popchange(3) = p.sigma_avg*pops(1);

    % cumulative P, cumulative A
    popchange(4) = p.rho*S;    
    popchange(5) = p.sigma_avg*pops(1); 

else

    % prescription classes
    popchange(1) = p.rho*S - (p.alpha*p.k + p.sigma_vec(1))*pops(1);
    for i = 2:p.k
        popchange(i) = p.alpha*p.k*pops(i-1) - (p.alpha*p.k + p.sigma_vec(i))*pops(i);
    end

    % susceptible, addicted
    popchange(p.k+1) = -p.rho*S + p.k*p.alpha*pops(p.k);
    popchange(p.k+2) =  sum(pops(1:p.k)'.*p.sigma_vec(1:end));

    % cumulative P, cumulative A
    popchange(p.k+3) = p.rho*S;
    popchange(p.k+4) = sum(pops(1:p.k)'.*p.sigma_vec(1:end));

end