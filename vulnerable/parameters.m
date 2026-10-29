function p = parameters(classes,range,sigavg);

p.range = range;
p.k = classes;

%% Model parameters

p.rho = 0.514;
p.alpha = 1/30*365;
p.sigma_avg = sigavg;
p.eta = 0.198;
p.delta = 0.7;
p.beta = 0.07;
p.beta_const = 0.00005;
p.fbase = 0.00744;

p.mu    = 0.0127;
p.muOD  = 0.0221;

p.gammaS2 = 12.0;
p.f0 = 0.00744;
p.fmax    = 0.02;
p.fmin    = 0.002;
p.cshift  = 365/12;
p.cstar   = 365/36;

if p.range > 0 && (p.k < 2 || p.range >= p.k*p.sigma_avg)
    disp('invalid range for prescription stages')
    return
end

%% Simulation setup

p.endT = 25;
p.N = 10000;

%% Prescription stages

p.num_classes = 2*p.k + 12;

p.sigma_vec = ones(p.k,1);

if p.range == 0
    p.sigma_vec = p.sigma_avg*p.sigma_vec';
else 
    fun = @(r) p.range * sum(r.^(0:p.k-1)) ...
        - p.k * p.sigma_avg * (r.^(p.k-1) - 1);
    upper = 2;
    while fun(upper) > 0
        upper = 2 * upper;
    end
    p.r = fzero(fun, [1 upper]);

    for q = 1:p.k
        temp = p.k.*p.sigma_avg.*p.r.^(q-1).*(1-p.r)/(1-p.r.^p.k);
        p.sigma_vec(q) = temp(1,1);
    end
    p.sigma_vec = p.sigma_vec'; 

    if isreal(p.r) == false
        disp('r is imaginary')
        p.r
        return
    elseif p.r <= 1
        disp('r is less than or equal to one')
        p.r
        return
    end

end

check_val = abs(mean(p.sigma_vec) - p.sigma_avg);
if check_val > 10e-9
    disp('average addiction rate does not match')
    p.k
    check_val
    return
end
