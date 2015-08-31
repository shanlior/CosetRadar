%function [stats] = sim2()
rng(1020);
g = global_settings();
N = 2;
for n=1:N
    targets = randomize_targets(g);    
%     targets.t = g.tau * [0.23 0.42 0.4201 0.72 0.84]';
%     targets.f = 1/g.tau * [0.875 0.69 0.2 0.105 0.33]';
    x = generate_analog_input_signal(g, targets);
    targets_MF_Nyq = classic_processing(g, x);
    targets_MF_SubNyq = classic_processing(g, x, g.sample_SubNyquist_factor, g.pulse_SubNyquist_factor);
    targets_CS_SCC = cs_processing_SCC2013(g, x);
    targets_CS_DF = cs_processing_DopplerFocusing(g, x);
    
    [targets_MF_Nyq, stats.MF_Nyq(n)] = analyze_results(g, targets, targets_MF_Nyq);
    [targets_MF_SubNyq, stats.MF_SubNyq(n)] = analyze_results(g, targets, targets_MF_SubNyq);
    [targets_CS_SCC, stats.CS_SCC(n)] = analyze_results(g, targets, targets_CS_SCC);
    [targets_CS_DF, stats.CS_DF(n)] = analyze_results(g, targets, targets_CS_DF);
    
    if mod(n,round(N/5)) == 0
        fprintf('n=%d, hit rate (MF-Nyq/MF-SubNyq/CS-SCC/CS-DF)=(%.2g,%.2g,%.2g,%.2g)\n', n, ...
            mean([stats.MF_Nyq.num_hits]) / g.L, ...
            mean([stats.MF_SubNyq.num_hits]) / g.L, ...
            mean([stats.CS_SCC.num_hits]) / g.L, ...
            mean([stats.CS_DF.num_hits]) / g.L);
    end    
%     plot_results(g, targets, targets_MF_Nyq, targets_MF_SubNyq, targets_CS_SCC, targets_CS_DF);
end

PrintStats('MF-Nyq   ', stats.MF_Nyq, g);
PrintStats('MF-SubNyq', stats.MF_SubNyq, g);
PrintStats('CS-SCC   ', stats.CS_SCC, g);
PrintStats('CS-DF    ', stats.CS_DF, g);
