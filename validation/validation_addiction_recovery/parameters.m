%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all model parameters %%
%% - All units are in years -  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = parameters(classes,range,sigavg);

p.range = range;
p.k = classes;

%% Model parameters

p.rho = 0.514; % battista 0.15   % cdc annual report, 2019
p.alpha = 1/30*365; % batistia 2.528  % rate return to to S, 18.4 days
p.sigma_avg = sigavg; % fixed addiction rate, doesn't distinguish between addiction types
p.eta = 0.198; % from lenhart

% Recovery rate


if p.k > 1 && p.range > p.k*p.sigma_avg
    disp('invalid range, range > k*sigma ')
    return
end


%% Simulation setup and demographic parameters

p.endT = 25;        % length of simulation (years)
p.N = 10000;%350000000;    % total population size, assumed constant


%% parameters for linear chain

if p.k == 1
    p.num_classes = p.k + 4;
else
    p.num_classes = 2*p.k + 4; % total number of classes (# prescribed classes, cumulative # addicted broken down by P class, S, A, cumulative S and A
end

p.sigma_vec = ones(p.k,1);

if p.k == 1
    p.sigma_vec = p.sigma_avg;
elseif p.range == 0
    p.sigma_vec = p.sigma_avg*p.sigma_vec';
else 
    myfun = @(x,t1,t2,t3) t2 * (1 - x.^t1) - t1 * t3 * (x.^(t1-1) - 1) * (1 - x);
    t1 = p.k;
    t2 = p.range;
    t3 = p.sigma_avg;

    fun = @(x) myfun(x,t1,t2,t3);
    x = fzero(fun,p.k);
    p.r = x;

    for q = 1:p.k
        temp = p.k.*p.sigma_avg.*p.r.^(q-1).*(1-p.r)/(1-p.r.^p.k);
        p.sigma_vec(q) = temp(1,1);
    end
    p.sigma_vec = p.sigma_vec'; 

    if isreal(p.r) == false
        disp('r is imaginary')
        p.r
        return
    elseif p.r < 1
        disp('r is less than one')
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

% a plot to look at: new addictions for each group over time. then look at
% sum. 



