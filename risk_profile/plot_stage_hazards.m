clear; close all;
k = 30;
sigma_avg = 0.00744;
multipliers = [0 1 2 5 10];

figure;
hold on;
for m = multipliers
    if m == 0
        sigma = sigma_avg * ones(1,k);
    else
        ratio = fzero(@(r) (r^(k-1)-1)/mean(r.^(0:k-1)) - m, [1+1e-8 2]);
        sigma = sigma_avg * ratio.^(0:k-1) / mean(ratio.^(0:k-1));
    end
    plot(1:k,sigma,'LineWidth',1.7,'DisplayName',sprintf('m_\\sigma = %g',m));
end
xlabel('Day of prescription');
ylabel('Addiction rate \sigma');
legend('Location','northwest');
xlim([1 k]);
