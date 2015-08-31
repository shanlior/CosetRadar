function [hit_rate, rms_t, bias_t, rms_f, bias_f] = sim3_collect_data(m, L, stats)

hit_rate = zeros(m, 7);
rms_t = zeros(m, 7);
bias_t = zeros(m, 7);
rms_f = zeros(m, 7);
bias_f = zeros(m, 7);
for k=1:m
    hit_rate(k, 1) = mean([stats(k).MF_Nyq.num_hits]) / L;
    hit_rate(k, 2) = mean([stats(k).MF_SubNyq.num_hits]) / L;
    hit_rate(k, 3) = mean([stats(k).CS_SCC.num_hits]) / L;
    hit_rate(k, 4) = mean([stats(k).DF_LPF.num_hits]) / L;
    hit_rate(k, 5) = mean([stats(k).DF_Random.num_hits]) / L;
    hit_rate(k, 6) = mean([stats(k).NU_SubNyq.num_hits]) / L;
    hit_rate(k, 7) = mean([stats(k).NU_SubNyq_Random.num_hits]) / L;
    
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).MF_Nyq, rms_t, bias_t, rms_f, bias_f, k, 1);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).MF_SubNyq, rms_t, bias_t, rms_f, bias_f, k, 2);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).CS_SCC, rms_t, bias_t, rms_f, bias_f, k, 3);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).DF_LPF, rms_t, bias_t, rms_f, bias_f, k, 4);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).DF_Random, rms_t, bias_t, rms_f, bias_f, k, 5);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).NU_SubNyq, rms_t, bias_t, rms_f, bias_f, k, 6);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).NU_SubNyq_Random, rms_t, bias_t, rms_f, bias_f, k, 7);
end
