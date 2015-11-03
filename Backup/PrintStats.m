function PrintStats(name, stats, g)

hit_rate = mean([stats.num_hits]) / g.L;
PDF = [stats.num_hits] / sum([stats.num_hits]);
rms_t = sqrt(mean_without_nans([stats.mse_t], PDF));
rms_f = sqrt(mean_without_nans([stats.mse_f], PDF));
bias_t = mean_without_nans([stats.bias_t], PDF);
bias_f = mean_without_nans([stats.bias_f], PDF);

fprintf('%s: hit_rate=%.2f RMS(t)=%.2f RMS(f)=%.2f bias(t)=%.2f bias(f)=%.2f\n', ...
    name, hit_rate, rms_t/g.Nyquist.dt, rms_f/g.Nyquist.df, bias_t/g.Nyquist.dt, bias_f/g.Nyquist.df);
