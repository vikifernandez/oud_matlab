clear; close all; clc;

%% Settings
maxP = 30;
sigavg = 0.00744;
range_vec = sigavg .* [0 1 2 5 10];
nc = length(range_vec);

%% Simulation and plot
figure('Position',[100 100 1000 600]);
for classy = 1:nc
    range = range_vec(classy);
    classes = maxP;
    p = parameters(classes,range,sigavg);
    p.beta = 0.0001;
    
    totalA  = zeros(maxP,1);
    prescA  = zeros(maxP,1);
    socialA_S   = zeros(maxP,1);
    socialA_P   = zeros(maxP,1);

    for j = 1:maxP
        p.cutoff = j;
        inits = zeros(p.num_classes,1); 
        prescribe_init = 0.0001;
        inits(1) = prescribe_init*p.N; 
        inits(p.k+2) = 0;
        inits(p.k+1) = p.N - sum(inits);
        
        options = odeset('NonNegative',1:p.num_classes); 
        [t, pop]=ode45(@opioid,[0 p.endT],inits,options,p);
        
        cumulA_total  = pop(:, p.k+4);
        cumulA_presc  = pop(:, 2*p.k + 5);
        cumulA_socialS = pop(:, 2*p.k + 6);
        cumulA_socialP = pop(:, 2*p.k + 7); 

        totalA(j)  = cumulA_total(end);
        prescA(j)  = cumulA_presc(end);
        socialA_S(j) = cumulA_socialS(end);
        socialA_P(j) = cumulA_socialP(end);
    end

    %% Plot risk profiles
    subplot(ceil(nc/2), 2, classy); hold on;

    area(1:maxP, [prescA, socialA_S, socialA_P], 'LineStyle','none');
    colormap([0.2 0.5 0.9; 0.98 0.65 0.25; 0.9 0.3 0.3]);
    alpha(0.55);

    plot(1:maxP, totalA, 'k', 'LineWidth', 1.6);

    xlabel('Prescription length (days)');
    ylabel('Cumulative addictions');
    title(sprintf('\\sigma-range = %.1fx mean \\sigma', range/sigavg));
    grid on;
    xlim([1 30]);

    if classy == 1
        legend({'Prescription-origin','Social-from-S','Social-from-P','Total'}, ...
               'Location','northwest');
    end
end

sgtitle('Cumulative addictions by source across σ-range values');
