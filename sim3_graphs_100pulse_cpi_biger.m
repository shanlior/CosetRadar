function sim3_graphs_100pulse_cpi_biger(g_LPF, Mp_all, hit_rate, rms_t, bias_t, rms_f, bias_f)

dir_name = fullfile('sim3_graphs_100pulse_cpi_biger','new2');
if ~exist(dir_name,'dir')
    mkdir(dir_name);
end

figure;
hold all;
% plot(Mp_all, hit_rate(:,1), '.-');
% plot(Mp_all, hit_rate(:,2), '*-');
% plot(Mp_all, hit_rate(:,3), 'o-');
% plot(Mp_all, hit_rate(:,4), 's-');
% plot(Mp_all, hit_rate(:,5), '^-');
plot(Mp_all, hit_rate(:,1), '.-');
plot(Mp_all, hit_rate(:,2), '*-');
title('Hit Rate');
add_stuff();
saveas(gcf, fullfile(dir_name, sprintf('fig1_%d.fig',length(Mp_all))));
saveas(gcf, fullfile(dir_name, sprintf('fig1_%d.jpg',length(Mp_all))));

figure;
subplot(211);
hold all;
% plot(Mp_all, rms_t(:,1) / g_LPF.Nyquist.dt, '.-');
% plot(Mp_all, rms_t(:,2) / g_LPF.Nyquist.dt, '*-');
% plot(Mp_all, rms_t(:,3) / g_LPF.Nyquist.dt, 'o-');
% plot(Mp_all, rms_t(:,4) / g_LPF.Nyquist.dt, 's-');
% plot(Mp_all, rms_t(:,5) / g_LPF.Nyquist.dt, '^-');
plot(Mp_all, rms_t(:,1) / g_LPF.Nyquist.dt, '.-');
plot(Mp_all, rms_t(:,2) / g_LPF.Nyquist.dt, '*-');
ylabel('[Nyquist bins]');
title('Time RMS error');
add_stuff();
subplot(212);
hold all;
% plot(Mp_all, rms_f(:,1) / g_LPF.Nyquist.df, '.-');
% plot(Mp_all, rms_f(:,2) / g_LPF.Nyquist.df, '*-');
% plot(Mp_all, rms_f(:,3) / g_LPF.Nyquist.df, 'o-');
% plot(Mp_all, rms_f(:,4) / g_LPF.Nyquist.df, 's-');
% plot(Mp_all, rms_f(:,5) / g_LPF.Nyquist.df, '^-');
plot(Mp_all, rms_f(:,1) / g_LPF.Nyquist.df, '.-');
plot(Mp_all, rms_f(:,2) / g_LPF.Nyquist.df, '*-');
ylabel('[Nyquist bins]');
title('Frequency RMS error');
add_stuff();
saveas(gcf, fullfile(dir_name, sprintf('fig2_%d.fig',length(Mp_all))));
saveas(gcf, fullfile(dir_name, sprintf('fig2_%d.jpg',length(Mp_all))));

figure;
subplot(211);
hold all;
% plot(Mp_all, bias_t(:,1) / g_LPF.Nyquist.dt, '.-');
% plot(Mp_all, bias_t(:,2) / g_LPF.Nyquist.dt, '*-');
% plot(Mp_all, bias_t(:,3) / g_LPF.Nyquist.dt, 'o-');
% plot(Mp_all, bias_t(:,4) / g_LPF.Nyquist.dt, 's-');
% plot(Mp_all, bias_t(:,5) / g_LPF.Nyquist.dt, '^-');
plot(Mp_all, bias_t(:,1) / g_LPF.Nyquist.dt, '.-');
plot(Mp_all, bias_t(:,2) / g_LPF.Nyquist.dt, '*-');
ylabel('[Nyquist bins]');
title('Time error bias');
add_stuff();
subplot(212);
hold all;
% plot(Mp_all, bias_f(:,1) / g_LPF.Nyquist.df, '.-');
% plot(Mp_all, bias_f(:,2) / g_LPF.Nyquist.df, '*-');
% plot(Mp_all, bias_f(:,3) / g_LPF.Nyquist.df, 'o-');
% plot(Mp_all, bias_f(:,4) / g_LPF.Nyquist.df, 's-');
% plot(Mp_all, bias_f(:,5) / g_LPF.Nyquist.df, '^-');
plot(Mp_all, bias_f(:,1) / g_LPF.Nyquist.df, '.-');
plot(Mp_all, bias_f(:,2) / g_LPF.Nyquist.df, '*-');
ylabel('[Nyquist bins]');
title('Frequency error bias');
add_stuff();
saveas(gcf, fullfile(dir_name, sprintf('fig3_%d.fig',length(Mp_all))));
saveas(gcf, fullfile(dir_name, sprintf('fig3_%d.jpg',length(Mp_all))));

function add_stuff()
% legend('MF-Nyquist','MF-SubNyquist','CS-LPF','CS-Random','CS-Blocks','CS-Cosets');
legend('Non-Uniform Sub-Nyquist - consecutive coeff. - FISTA', 'Non-Uniform Sub-Nyquist - random coeff. - FISTA');
xlabel('CPI');
grid on;
