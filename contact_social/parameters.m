function p = parameters(classes,range,sigavg);

p.range = range;
p.k = classes;

%% Model parameters

p.rho = 0.514;
p.alpha = 1/30*365;
p.sigma_avg = sigavg;
p.eta = 0.198;

p.beta = 0.07;
p.mu    = 0.0127;
p.muOD  = 0.0221;

if p.k > 1 && p.range > p.k*p.sigma_avg
    disp('invalid range, range > k*sigma ')
    return
end

%% Simulation setup

p.endT = 25;
p.N = 10000;

%% Prescription stages

p.num_classes = 2*p.k + 7;

p.sigma_vec = ones(p.k,1);

if p.range == 0
    p.sigma_vec = p.sigma_avg*p.sigma_vec';
else 
    myfun = @(x,t1,t2,t3) t2 * (1 - x.^t1) - t1 * t3 * (x.^(t1-1) - 1) * (1 - x);
    t1 = p.k;
    t2 = p.range;
    t3 = p.sigma_avg;

    fun = @(x) myfun(x,t1,t2,t3);
    x = fsolve(fun,2);
    p.r = x;
    if x == NaN
        x = fzero(fun,1.05);
        p.r = x;
    end
    p.r

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
