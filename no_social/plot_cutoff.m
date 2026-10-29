clear; close all; clc;

%% Settings
maxP = 30;
sigavg = 0.00744;
range_vec = sigavg .* [0 1 2 5 10];
nc = numel(range_vec);

%% Risk profiles
figure('Position', [100 200 700 900]);
hold on;
for classy = 1:nc

    classes = maxP;
    range = range_vec(classy);
    p = parameters(classes,range,sigavg);
    totalA = zeros(maxP,1);

    for j = 1:maxP

        p.cutoff = j;
        inits = zeros(p.num_classes,1);
        prescribe_init = 0.0001;
        inits(1) = prescribe_init*p.N;
        inits(p.k+2) = 0;
        inits(p.k+1) = p.N - sum(inits);

       options = odeset('NonNegative',1:p.num_classes); 
       [~, pop]=ode45(@opioid,[0 p.endT],inits,options,p);

        cumulA = pop(:,p.k+4);

        totalA(j) = cumulA(end);
     
    end

    plot(1:maxP, totalA, 'LineWidth', 1.4);

end

xlabel('Prescription length (days)');
ylabel('Cumulative addictions');
legend({'m_{\sigma} = 0', 'm_{\sigma} = 1', 'm_{\sigma} = 2', ...
        'm_{\sigma} = 5', 'm_{\sigma} = 10'}, 'Location','northwest');
