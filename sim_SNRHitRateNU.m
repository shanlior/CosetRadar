function [stats] = sim_SNRHitRateFISTA() 
rng(0);
SNR_dB_all = -40:1:-39; % for P=100
N = 1;
L = 5;
P = 100;
CoeffDist = 0; % CoeffDist=0 - consecutive coeff; CoeffDist=1 - random coeff
datestr(now)
tic;

stats.NU_100 = repmat(create_stats(), N, 1);
stats.NU_50  = repmat(create_stats(), N, 1);
stats.NU_20  = repmat(create_stats(), N, 1);
stats.NU_10  = repmat(create_stats(), N, 1);
stats = repmat(stats, length(SNR_dB_all), 1);

for m=1:length(SNR_dB_all)
    disp(['**** SNR = ' num2str(SNR_dB_all(m)) '[dB] ****']);
   
    for n=1:N
        g_NU_100 = global_settings_nu(100 ,P, L, SNR_dB_all(m), 10, 1, CoeffDist, 1);
        g_NU_50  = global_settings_nu(50,  P, L, SNR_dB_all(m), 10, 1, CoeffDist, 1);
        g_NU_20  = global_settings_nu(20,  P, L, SNR_dB_all(m), 10, 1, CoeffDist, 1);
        g_NU_10  = global_settings_nu(10,  P, L, SNR_dB_all(m), 10, 1, CoeffDist, 1);
        
        targets = randomize_targets(g_NU_100);  % Doesn't matter which global settings
        x = generate_analog_input_signal(g_NU_100, targets);  
        targets_NU_100 = non_uniform_sub_nyquist(g_NU_100, x, 1);
        targets_NU_50  = non_uniform_sub_nyquist(g_NU_50, x,  1);
        targets_NU_20  = non_uniform_sub_nyquist(g_NU_20, x,  1);
        targets_NU_10  = non_uniform_sub_nyquist(g_NU_10, x,  1);
        
        [targets_NU_100, stats(m).NU_100(n)] = analyze_results(g_NU_100, targets, targets_NU_100);
        [targets_NU_50,  stats(m).NU_50(n)]  = analyze_results(g_NU_50,  targets, targets_NU_50);
        [targets_NU_20,  stats(m).NU_20(n)]  = analyze_results(g_NU_20,  targets, targets_NU_20);
        [targets_NU_10,  stats(m).NU_10(n)]  = analyze_results(g_NU_10,  targets, targets_NU_10);
        
        if mod(n,round(N/5)) == 0
            fprintf('n=%d, hit rate (MF-Nyq/MF-SubNyq/CS-SCC/CS-DF-LPF/CS-DF-Random)=(%.2g,%.2g,%.2g,%.2g,%.2g)\n', n, ...
                mean([stats(m).NU_100(1:n).num_hits]) / g_NU_100.L, ...
                mean([stats(m).NU_50(1:n).num_hits])  / g_NU_50.L, ...
                mean([stats(m).NU_20(1:n).num_hits])  / g_NU_20.L, ...
                mean([stats(m).NU_10(1:n).num_hits])  / g_NU_10.L);
        end
    end
    PrintStats('NU-100/100', stats(m).NU_100, g_NU_100);
    PrintStats('NU-50/100',  stats(m).NU_50,  g_NU_50);
    PrintStats('NU-20/100',  stats(m).NU_20,  g_NU_20);
    PrintStats('NU-10/100',  stats(m).NU_10,  g_NU_10);
    
    [hit_rate, rms_t, bias_t, rms_f, bias_f] = sim_SNRHitRateNU_collect_data(m, L, stats);
    close all;
    sim_SNRHitRateNU_graphs(global_settings(), SNR_dB_all(1:m), hit_rate, rms_t, bias_t, rms_f, bias_f);
end

toc;
datestr(now)
if CoeffDist==0
    save('SNRHitRateFistaConCoeff');
else
    save('SNRHitRateFistaRandCoeff');
end
