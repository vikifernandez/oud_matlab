clear; close all;
k = 30;
sigma_avg = 0.00744;
multipliers = [0 1 2 5 10];
basic = zeros(k,numel(multipliers));
first_time = basic;
relapses = basic;

for q = 1:numel(multipliers)
    p = parameters(k,multipliers(q)*sigma_avg,sigma_avg);
    p.num_classes = 2*k+10;
    for cutoff = 1:k
        p.cutoff = cutoff;
        p.alpha = 365/p.k;
        initial = zeros(p.num_classes,1);
        initial(1) = 0.0001*p.N;
        initial(k+1) = p.N-initial(1);
        options = odeset('NonNegative',1:p.num_classes);
        [~,result] = ode45(@opioid_basic,[0 p.endT],initial,options,p);
        basic(cutoff,q) = result(end,k+4);
        [~,result] = ode45(@opioid_relapse,[0 p.endT],initial,options,p);
        first_time(cutoff,q) = result(end,2*k+5);
        relapses(cutoff,q) = result(end,2*k+6);
    end
end

figure('Position',[100 100 1600 900]);
subplot(2,2,1); hold on; grid on;
for q = 1:numel(multipliers)
    p = parameters(k,multipliers(q)*sigma_avg,sigma_avg);
    plot(1:k,p.sigma_vec,'LineWidth',1.5, ...
        'DisplayName',sprintf('\\sigma-range = %.1fx',multipliers(q)));
end
title('Addiction hazard vs prescription day');
xlabel('Prescription Day'); ylabel('\sigma_i'); legend('Location','northwest');

titles = {'Basic model','Relapse model: first-time addictions', ...
    'Relapse model: first-time (solid) vs relapse (dotted)'};
for panel = 2:4
    subplot(2,2,panel); hold on; grid on;
    for q = 1:numel(multipliers)
        color = get(gca,'ColorOrder'); color = color(q,:);
        if panel == 2, y = basic(:,q); else, y = first_time(:,q); end
        plot(1:k,y,'Color',color,'LineWidth',1.5);
        if panel == 4
            plot(1:k,relapses(:,q),'--','Color',color,'LineWidth',1.5);
        end
    end
    title(titles{panel-1});
    xlabel('Prescription cutoff (days)');
    ylabel('Cumulative addictions');
end
