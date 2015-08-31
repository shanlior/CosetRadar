function [rms_t, bias_t, rms_f, bias_f] = calc_rms_bias(stat, rms_t, bias_t, rms_f, bias_f, m, index)
PDF = [stat.num_hits] / sum([stat.num_hits]);
rms_t(m, index) = sqrt(mean_without_nans([stat.mse_t], PDF));
bias_t(m, index) = mean_without_nans([stat.bias_t], PDF);
rms_f(m, index) = sqrt(mean_without_nans([stat.mse_f], PDF));
bias_f(m, index) = mean_without_nans([stat.bias_f], PDF);