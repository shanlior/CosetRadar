function [stats] = sim3()
rng(0);
% SNR_dB_all = -30:1:10; % for P=1
% SNR_dB_all = -50:2:0; % for P=100
SNR_dB_all = -15:1:-10; % for P=100
N = 2;
L = 5;
P = 100;
datestr(now)
tic;

stats.MF_Nyq = repmat(create_stats(), N, 1);
stats.MF_SubNyq = repmat(create_stats(), N, 1);
stats.CS_SCC = repmat(create_stats(), N, 1);
stats.DF_LPF = repmat(create_stats(), N, 1);
stats.DF_Random = repmat(create_stats(), N, 1);
stats.NU_SubNyq = repmat(create_stats(), N, 1);
stats.NU_SubNyq_Random = repmat(create_stats(), N, 1);
stats = repmat(stats, length(SNR_dB_all), 1);

for m=1:length(SNR_dB_all)
    disp(['**** SNR = ' num2str(SNR_dB_all(m)) '[dB] ****']);
    g_LPF = global_settings(P, L, SNR_dB_all(m), 10, 1, 0, 1);
    g_Random = global_settings(P, L, SNR_dB_all(m), 10, 1, 1, 1);
    g_LPF_nu = global_settings_nu(50 ,P, L, SNR_dB_all(m), 10, 1, 0, 1);
    g_Random_nu = global_settings_nu(50, P, L, SNR_dB_all(m), 10, 1, 1, 1);
    for n=1:N
        targets = randomize_targets(g_LPF);
        x = generate_analog_input_signal(g_LPF, targets);
        targets_MF_Nyq = classic_processing(g_LPF, x);
        targets_MF_SubNyq = classic_processing(g_LPF, x, g_LPF.sample_SubNyquist_factor, g_LPF.pulse_SubNyquist_factor);
        targets_CS_SCC = cs_processing_SCC2013(g_LPF, x);
        targets_DF_LPF = cs_processing_DopplerFocusing(g_LPF, x);
        targets_DF_Random = cs_processing_DopplerFocusing(g_Random, x);
        targets_NU_SubNyq = non_uniform_sub_nyquist(g_LPF_nu, x);
        targets_NU_SubNyq_Random = non_uniform_sub_nyquist(g_Random_nu, x);
        
%         plot_results(g, targets, targets_MF_Nyq, targets_MF_SubNyq, targets_CS_SCC, targets_CS_DF);
        
        [targets_MF_Nyq, stats(m).MF_Nyq(n)] = analyze_results(g_LPF, targets, targets_MF_Nyq);
        [targets_MF_SubNyq, stats(m).MF_SubNyq(n)] = analyze_results(g_LPF, targets, targets_MF_SubNyq);
        [targets_CS_SCC, stats(m).CS_SCC(n)] = analyze_results(g_LPF, targets, targets_CS_SCC);
        [targets_DF_LPF, stats(m).DF_LPF(n)] = analyze_results(g_LPF, targets, targets_DF_LPF);
        [targets_DF_Random, stats(m).DF_Random(n)] = analyze_results(g_Random, targets, targets_DF_Random);
        [targets_NU_SubNyq, stats(m).NU_SubNyq(n)] = analyze_results(g_LPF, targets, targets_NU_SubNyq);
        [targets_NU_SubNyq_Random, stats(m).NU_SubNyq_Random(n)] = analyze_results(g_Random, targets, targets_NU_SubNyq_Random);
        
        if mod(n,round(N/5)) == 0
            fprintf('n=%d, hit rate (MF-Nyq/MF-SubNyq/CS-SCC/CS-DF-LPF/CS-DF-Random)=(%.2g,%.2g,%.2g,%.2g,%.2g)\n', n, ...
                mean([stats(m).MF_Nyq(1:n).num_hits]) / g_LPF.L, ...
                mean([stats(m).MF_SubNyq(1:n).num_hits]) / g_LPF.L, ...
                mean([stats(m).CS_SCC(1:n).num_hits]) / g_LPF.L, ...
                mean([stats(m).DF_LPF(1:n).num_hits]) / g_LPF.L, ...
                mean([stats(m).DF_Random(1:n).num_hits]) / g_Random.L, ...
                mean([stats(m).NU_SubNyq(1:n).num_hits]) / g_LPF.L, ...
                mean([stats(m).NU_SubNyq_Random(1:n).num_hits]) / g_Random.L);
        end
    end
    PrintStats('MF-Nyq   ', stats(m).MF_Nyq, g_LPF);
    PrintStats('MF-SubNyq', stats(m).MF_SubNyq, g_LPF);
    PrintStats('CS-SCC   ', stats(m).CS_SCC, g_LPF);
    PrintStats('DF-LPF   ', stats(m).DF_LPF, g_LPF);
    PrintStats('DF-Random', stats(m).DF_Random, g_Random);
    PrintStats('NU-SubNyq', stats(m).NU_SubNyq, g_LPF);
    PrintStats('NU-SubNyq-Random', stats(m).NU_SubNyq_Random, g_Random);
    
    [hit_rate, rms_t, bias_t, rms_f, bias_f] = sim3_collect_data(m, L, stats);
    close all;
    sim3_graphs(global_settings(), SNR_dB_all(1:m), hit_rate, rms_t, bias_t, rms_f, bias_f);
end
toc;
datestr(now)
save;
