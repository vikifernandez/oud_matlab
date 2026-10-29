clear; close all; clc;

%% Setup

maxP   = 30;
sigavg = 0.00744;
mults  = [0 1 2 5 10];
nm     = length(mults);

cutoff_vec = 1:maxP;
cvals = cutoff_vec / 365;

u_mid_day   = (cutoff_vec(1) + cutoff_vec(end)) / 2;
u_quad_vals = (cutoff_vec - u_mid_day).^2;
u_quad_bar  = mean(u_quad_vals);

forms = {'increasing','decreasing','ushaped'};
form_labels = {'Increasing f(c)', 'Decreasing f(c)', 'U-shaped f(c)'};
nf = length(forms);

spread_mults = [0.25 0.5 0.75 1 1.25 1.5 2];
ns = length(spread_mults);

base_Adec = 0.006;
base_Ainc = 0.006;

base_mU   = 1.82e-5;

opt_cutoff = zeros(ns, nm, nf);

%% Simulation

for mm = 1:nm

    fprintf('\nRunning sigma multiplier %.2f...\n', mults(mm));

    range = mults(mm) * sigavg;

    p = parameters(maxP, range, sigavg);

    %% Post-prescription parameters
    p.gammaS2 = 12;
    p.f0      = 0.00744;
    p.fbase   = 0.00744;

    p.cshift = 12/365;

    p.u_mid_day  = u_mid_day;
    p.u_quad_bar = u_quad_bar;

    gdec_vals = p.cshift ./ (cvals + p.cshift);
    ginc_vals = cvals    ./ (cvals + p.cshift);

    p.gdec_bar = mean(gdec_vals);
    p.ginc_bar = mean(ginc_vals);

    for ss = 1:ns

        spread = spread_mults(ss);

        fprintf('  Spread %.2f...\n', spread);

        p.Adec = spread * base_Adec;
        p.Ainc = spread * base_Ainc;
        p.mU   = spread * base_mU;

        fdec_check = p.fbase + p.Adec * (gdec_vals - p.gdec_bar);
        finc_check = p.fbase + p.Ainc * (ginc_vals - p.ginc_bar);
        fU_check   = p.fbase + p.mU * (u_quad_vals - p.u_quad_bar);

        if min([fdec_check finc_check fU_check]) < 0
            warning('Some f(c) values are negative for sigma %.2f, spread %.2f.', ...
                mults(mm), spread);
        end

        for ff = 1:nf

            p.f_form = forms{ff};

            total_burden = zeros(maxP,1);

            %% Prescription cutoffs
            for j = 1:maxP

                p.cutoff = j;
                p.c      = j/365;

                p.num_classes = p.k + 11;

                inits = zeros(p.num_classes,1);
                prescribe_init = 0.0001 * p.N;

                inits(1)     = prescribe_init;
                inits(p.k+1) = p.N - prescribe_init;
                inits(p.k+2) = 0;
                inits(p.k+3) = 0;
                inits(p.k+4) = 0;

                options = odeset('NonNegative',1:p.num_classes);

                [~,pop] = ode45(@NEWu_opioid_postvuln_social, [0 p.endT], inits, options, p);

                first_time = pop(end,p.k+6);
                relapse    = pop(end,p.k+7);

                total_burden(j) = first_time + relapse;
            end

            [~, idx_min] = min(total_burden);
            opt_cutoff(ss,mm,ff) = cutoff_vec(idx_min);

        end
    end
end

%% Optimal cutoffs and f(c)
figure('Position',[100 100 1700 500]);
for ff = 1:nf
    subplot(1,nf,ff);
    imagesc(mults,spread_mults,opt_cutoff(:,:,ff));
    set(gca,'YDir','normal');
    xlabel('m_{\sigma} multiplier');
    ylabel('f(c) spread multiplier');
    title(form_labels{ff});
    cb = colorbar;
    cb.Label.String = 'Optimal cutoff c^* (days)';
    caxis([1 maxP]);
    grid off;
end
sgtitle('Optimal prescription cutoff as \\sigma-range and f(c) spread vary');

figure('Position',[100 100 1400 400]);
colors = parula(numel(spread_mults)+1);
for ff = 1:nf
    subplot(1,nf,ff); hold on; grid on; box on;
    for ss = 1:ns
        spread = spread_mults(ss);
        switch forms{ff}
            case 'decreasing'
                fvals = sigavg + spread*base_Adec*(gdec_vals-mean(gdec_vals));
            case 'increasing'
                fvals = sigavg + spread*base_Ainc*(ginc_vals-mean(ginc_vals));
            case 'ushaped'
                fvals = sigavg + spread*base_mU*(u_quad_vals-u_quad_bar);
        end
        plot(cutoff_vec,fvals,'Color',colors(ss+1,:), ...
            'LineWidth',2,'DisplayName',sprintf('s = %.2f',spread));
    end
    yline(sigavg,'--k','LineWidth',1.5,'DisplayName','fbase');
    xlabel('Prescription cutoff (days)');
    ylabel('f(c)');
    title(form_labels{ff});
    xlim([1 maxP]); ylim([0 0.017]);
    xticks([1 5 10 15 20 25 30]);
    yticks([0 0.005 0.010 0.015]);
    if ff == 1, legend('Location','best'); end
end
sgtitle('Diagnostic: f(c) curves remain centered at fbase = 0.00744');
