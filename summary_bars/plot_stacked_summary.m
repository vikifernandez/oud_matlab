function plot_stacked_summary
close all;

model_dir = fullfile(fileparts(mfilename('fullpath')),'..','combined_social');
addpath(model_dir,'-begin');
path_cleanup = onCleanup(@() rmpath(model_dir));

k = 30;
sigma_avg = 0.00744;
labels = {'Constant','Contact-dependent','Both'};
constant_rates = [0.0001 0 0.00005];
contact_rates = [0 0.07 0.07];
basic = zeros(1,3);
first_time = zeros(1,3);
relapses = zeros(1,3);

for variant = 1:3
    p = parameters(k,10*sigma_avg,sigma_avg);
    p.cutoff = 30;
    p.alpha = 365/30;
    p.beta_const = constant_rates(variant);
    p.beta = contact_rates(variant);
    p.num_classes = 2*k+10;

    initial = zeros(p.num_classes,1);
    initial(1) = 0.0001*p.N;
    initial(k+1) = p.N-initial(1);
    options = odeset('NonNegative',1:p.num_classes);

    [~,states] = ode45(@opioid_basic,[0 p.endT],initial,options,p);
    basic(variant) = states(end,k+4);
    [~,states] = ode45(@opioid_relapse,[0 p.endT],initial,options,p);
    first_time(variant) = states(end,2*k+5);
    relapses(variant) = states(end,2*k+6);
end

x = 1:3;
figure('Units','pixels','Position',[100 100 1000 900]);
h_basic = bar(x-0.18,basic,0.32,'FaceColor',[0 0.447 0.741]);
hold on;
h_relapse = bar(x+0.18,[first_time(:) relapses(:)],0.32,'stacked');
h_relapse(1).FaceColor = [0.85 0.325 0.098];
h_relapse(2).FaceColor = [0.929 0.694 0.125];
set(gca,'XTick',x,'XTickLabel',labels,'XTickLabelRotation',15);
ylabel('Cumulative addictions');
legend([h_basic h_relapse(1) h_relapse(2)], ...
    {'Basic model','Relapse model: new addictions','Relapse model: relapses'}, ...
    'Location','northwest');
grid on;
set(gca, 'FontSize', 14);  
ylabel('Cumulative addictions', 'FontSize', 14);
legend('FontSize', 14);   
end
