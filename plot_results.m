function plot_results(g, targets, targets_NU_SubNyq)

figure;
cla;
if g.hit_rate_threshold.f > 1/g.tau
    xlabel('time [sec]');
    axis([0 g.Q*g.tau -1 1]);
    grid on;
    hold on;
    for l=1:g.L % draw hit rate ellipses
        h_real=plot(targets.t(l), 0, 'b.');
        plot(targets.t(l) + g.hit_rate_threshold.t * [-1 1 1 -1 -1], 0.1*[-1 -1 1 1 -1], 'b');
    end
    h_NU_SubNyquist = Draw1D(targets, targets_NU_SubNyq, 's', 'cyan'); 
else
    xlabel('time [sec]');
    ylabel('frequency [hz]');
    axis([0 (g.Q*g.tau) 0 1/g.tau]);
    grid on;
    hold on;
    theta = linspace(0, 2*pi, 40);
    for l=1:g.L % draw hit rate ellipses
        h_real=plot(targets.t(l), targets.f(l), 'b.');
        for k1=-1:1
            for k2=-1:1
                plot(k1*g.tau + targets.t(l) + g.hit_rate_threshold.t * cos(theta), ...
                    k2/g.tau + targets.f(l) + g.hit_rate_threshold.f * sin(theta), 'b');
            end
        end
    end
    
    h_NU_SubNyquist = Draw2D(targets, targets_NU_SubNyq, 's', 'cyan');
end
legend([h_NU_SubNyquist(1)], ...
    { 'CosetRadar'}, 'Location', 'Best');
snr_db = 10*log10(g.snr);
title_str{1} = sprintf('L=%d P=%d Avg. Sample SNR=%.1f[dB]', g.L, g.P, snr_db);
title_str{2} = sprintf('UnderSampling ratio sample/pulse=(%.1f,%.1f)', g.sample_SubNyquist_factor, ...
    g.pulse_SubNyquist_factor);
% title_str{3} = sprintf(' %.2f (NU-SubNyq)', ...
%      mean(targets_NU_SubNyq.best_perm_hits));
title(title_str);
hold off;

function [h] = Draw1D(targets, targets_est, center_marker, color)

h = [];
for l=1:length(targets.a)
    est_t = targets_est.t(targets_est.best_perm(l));
    h = [h plot(est_t, 0, center_marker)];
    if targets_est.best_perm_hits(l)
        h1 = line([targets.t(l) est_t],[0 0]);
        set(h1,'Color',color,'LineStyle',':');
    end
end

function [h] = Draw2D(targets, targets_est, center_marker, color)

h = [];
for l=1:length(targets.a)
    est_t = targets_est.t(targets_est.best_perm(l));
    est_f = targets_est.f(targets_est.best_perm(l));
%     est_t = targets_est.t(l);
%     est_f = targets_est.f(l);
    h = [h plot(est_t, est_f, center_marker)];
    if targets_est.best_perm_hits(l)
%     if 0
        h1 = line([targets.t(l) est_t],[targets.f(l) est_f]);
        set(h1,'Color',color,'LineStyle',':');
    end
end
