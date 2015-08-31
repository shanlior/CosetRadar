function sim3_graphs(g_LPF, SNR_dB_all, hit_rate, rms_t, bias_t, rms_f, bias_f)

dir_name = fullfile('sim3_graphs','new1');
if ~exist(dir_name,'dir')
    mkdir(dir_name);
end

figure;
hold all;
plot(SNR_dB_all, hit_rate(:,1), '.-');
plot(SNR_dB_all, hit_rate(:,2), '*-');
plot(SNR_dB_all, hit_rate(:,3), 'o-');
plot(SNR_dB_all, hit_rate(:,4), 's-');
plot(SNR_dB_all, hit_rate(:,5), '^-');
plot(SNR_dB_all, hit_rate(:,6), '>-');
plot(SNR_dB_all, hit_rate(:,7), '<-');
title('Hit Rate');
add_stuff();
saveas(gcf, fullfile(dir_name, sprintf('fig1_%d.fig',length(SNR_dB_all))));
saveas(gcf, fullfile(dir_name, sprintf('fig1_%d.jpg',length(SNR_dB_all))));

figure;
subplot(211);
hold all;
plot(SNR_dB_all, rms_t(:,1) / g_LPF.Nyquist.dt, '.-');
plot(SNR_dB_all, rms_t(:,2) / g_LPF.Nyquist.dt, '*-');
plot(SNR_dB_all, rms_t(:,3) / g_LPF.Nyquist.dt, 'o-');
plot(SNR_dB_all, rms_t(:,4) / g_LPF.Nyquist.dt, 's-');
plot(SNR_dB_all, rms_t(:,5) / g_LPF.Nyquist.dt, '^-');
plot(SNR_dB_all, rms_t(:,6) / g_LPF.Nyquist.dt, '>-');
plot(SNR_dB_all, rms_t(:,7) / g_LPF.Nyquist.dt, '<-');
ylabel('[Nyquist bins]');
title('Time RMS error');
add_stuff();
subplot(212);
hold all;
plot(SNR_dB_all, rms_f(:,1) / g_LPF.Nyquist.df, '.-');
plot(SNR_dB_all, rms_f(:,2) / g_LPF.Nyquist.df, '*-');
plot(SNR_dB_all, rms_f(:,3) / g_LPF.Nyquist.df, 'o-');
plot(SNR_dB_all, rms_f(:,4) / g_LPF.Nyquist.df, 's-');
plot(SNR_dB_all, rms_f(:,5) / g_LPF.Nyquist.df, '^-');
plot(SNR_dB_all, rms_f(:,6) / g_LPF.Nyquist.df, '>-');
plot(SNR_dB_all, rms_f(:,7) / g_LPF.Nyquist.df, '<-');
ylabel('[Nyquist bins]');
title('Frequency RMS error');
add_stuff();
saveas(gcf, fullfile(dir_name, sprintf('fig2_%d.fig',length(SNR_dB_all))));
saveas(gcf, fullfile(dir_name, sprintf('fig2_%d.jpg',length(SNR_dB_all))));

figure;
subplot(211);
hold all;
plot(SNR_dB_all, bias_t(:,1) / g_LPF.Nyquist.dt, '.-');
plot(SNR_dB_all, bias_t(:,2) / g_LPF.Nyquist.dt, '*-');
plot(SNR_dB_all, bias_t(:,3) / g_LPF.Nyquist.dt, 'o-');
plot(SNR_dB_all, bias_t(:,4) / g_LPF.Nyquist.dt, 's-');
plot(SNR_dB_all, bias_t(:,5) / g_LPF.Nyquist.dt, '^-');
plot(SNR_dB_all, bias_t(:,6) / g_LPF.Nyquist.dt, '>-');
plot(SNR_dB_all, bias_t(:,7) / g_LPF.Nyquist.dt, '<-');
ylabel('[Nyquist bins]');
title('Time error bias');
add_stuff();
subplot(212);
hold all;
plot(SNR_dB_all, bias_f(:,1) / g_LPF.Nyquist.df, '.-');
plot(SNR_dB_all, bias_f(:,2) / g_LPF.Nyquist.df, '*-');
plot(SNR_dB_all, bias_f(:,3) / g_LPF.Nyquist.df, 'o-');
plot(SNR_dB_all, bias_f(:,4) / g_LPF.Nyquist.df, 's-');
plot(SNR_dB_all, bias_f(:,5) / g_LPF.Nyquist.df, '^-');
plot(SNR_dB_all, bias_f(:,6) / g_LPF.Nyquist.df, '>-');
plot(SNR_dB_all, bias_f(:,7) / g_LPF.Nyquist.df, '<-');
ylabel('[Nyquist bins]');
title('Frequency error bias');
add_stuff();
saveas(gcf, fullfile(dir_name, sprintf('fig3_%d.fig',length(SNR_dB_all))));
saveas(gcf, fullfile(dir_name, sprintf('fig3_%d.jpg',length(SNR_dB_all))));

function add_stuff()
% legend('MF-Nyquist','MF-SubNyquist','CS-LPF','CS-Random','CS-Blocks','CS-Cosets');
legend('Classic processing(Nyquist)', 'Classic processing(Sub-Nyquist)', ...
    'Two stage CS(Sub-Nyq)', 'Doppler Focusing - consecutive coeff.(Sub-Nyq)', ...
    'Doppler Focusing - random coeff.(Sub-Nyq)', ...
    'Non-Uniform Sub-Nyquist - consecutive coeff. - FISTA', 'Non-Uniform Sub-Nyquist - random coeff. - FISTA');
xlabel('SNR [dB]');
grid on;
