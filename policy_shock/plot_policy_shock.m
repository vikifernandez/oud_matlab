clear; close all; clc;

%% Setup
maxP   = 15;
sigavg = 0.00744;

multipliers = [0 1 2 5];
t1 = linspace(0, 1, 60);
t2 = linspace(1, 25, 690);

DEN_THRESH = 1e-6;

figure('Position',[100 100 1700 450]);

for i = 1:numel(multipliers)

    mult  = multipliers(i);
    range = sigavg * mult;

    scenarios = [3 7];
    out = struct();

    for s = 1:2
        cutoff2 = scenarios(s);

        %% Phase 1
        p = parameters15(maxP, range, sigavg);
        p.cutoff = 15;
        p.alpha  = 1/15 * 365;

        p.num_classes = 2*p.k + 11;

        inits = zeros(p.num_classes,1);
        prescribe_init = 0.0001;

        inits(1)      = prescribe_init * p.N;
        inits(p.k+2)  = 0;
        inits(p.k+3)  = 0;
        inits(p.k+1)  = p.N - sum(inits);

        options = odeset('NonNegative',1:p.num_classes);

        [t_phase1, pop1] = ode45(@opioid_relapse, t1, inits, options, p);

        %% Phase 2
        p.cutoff = cutoff2;

        inits2 = pop1(end,:)';
        [t_phase2, pop2] = ode45(@opioid_relapse, t2, inits2, options, p);

        %% Combine
        t_full   = [t_phase1 ; t_phase2];
        pop_full = [pop1 ; pop2];

        %% Cumulative series
        C_first   = pop_full(:, 2*p.k + 5);
        C_relapse = pop_full(:, 2*p.k + 6);
        C_add_tot = C_first + C_relapse;
        C_OD      = pop_full(:, 2*p.k + 11);

        out(s).cutoff = cutoff2;
        out(s).t      = t_full;
        out(s).C_add  = C_add_tot;
        out(s).C_first= C_first;
        out(s).C_OD   = C_OD;
    end

    idx3 = find([out.cutoff] == 3, 1);
    idx7 = find([out.cutoff] == 7, 1);

    t = out(idx7).t;

    base_add   = out(idx7).C_add;
    base_first = out(idx7).C_first;
    base_od    = out(idx7).C_OD;

    pct_add   = 100 * (out(idx3).C_add   - base_add)   ./ base_add;
    pct_first = 100 * (out(idx3).C_first - base_first) ./ base_first;
    pct_od    = 100 * (out(idx3).C_OD    - base_od)    ./ base_od;

    pct_add(base_add   < DEN_THRESH) = NaN;
    pct_first(base_first < DEN_THRESH) = NaN;
    pct_od(base_od     < DEN_THRESH) = NaN;

    subplot(1,4,i); hold on; grid on;

    h = plot(t, pct_add, 'LineWidth', 2.2, ...
        'DisplayName','%Δ Cumulative addictions (total)');

    plot(t, pct_first, '--', 'Color', h.Color, 'LineWidth', 2.2, ...
        'DisplayName','%Δ First-time addictions');

    plot(t, pct_od, ':', 'Color', h.Color, 'LineWidth', 2.4, ...
        'DisplayName','%Δ Cumulative overdose deaths');

    xline(1,'k:','Policy change', ...
    'LineWidth',1.5, ...
    'LabelVerticalAlignment','middle', ...
    'LabelHorizontalAlignment','center', ...
    'HandleVisibility','off');
    yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');

    title(sprintf('%%Δ (3 vs 7), multiplier = %gx', mult));
    xlabel('Time (years)');
    if i == 1
        ylabel('% difference vs 7-day cutoff');
    end

    if i == 4
        legend('Location','northwest');
    end
end

axs = findall(gcf,'Type','axes');
Ymax = -inf; Ymin = inf;
for a = axs'
    yl = a.YLim;
    Ymin = min(Ymin, yl(1));
    Ymax = max(Ymax, yl(2));
end
for a = axs'
    ylim(a, [Ymin Ymax]);
end
