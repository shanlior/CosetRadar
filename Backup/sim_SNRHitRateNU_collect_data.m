function [hit_rate, rms_t, bias_t, rms_f, bias_f] = sim_SNRHitRateFISTA_collect_data(m, L, stats)

hit_rate = zeros(m, 4);
rms_t = zeros(m, 4);
bias_t = zeros(m, 4);
rms_f = zeros(m, 4);
bias_f = zeros(m, 4);

for k=1:m
    hit_rate(k, 1) = mean([stats(k).NU_100.num_hits]) / L;
    hit_rate(k, 2) = mean([stats(k).NU_50.num_hits])  / L;
    hit_rate(k, 3) = mean([stats(k).NU_20.num_hits])  / L;
    hit_rate(k, 4) = mean([stats(k).NU_10.num_hits])  / L;

    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).NU_100, rms_t, bias_t, rms_f, bias_f, k, 1);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).NU_50,  rms_t, bias_t, rms_f, bias_f, k, 2);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).NU_20,  rms_t, bias_t, rms_f, bias_f, k, 3);
    [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stats(k).NU_10,  rms_t, bias_t, rms_f, bias_f, k, 4);
end
