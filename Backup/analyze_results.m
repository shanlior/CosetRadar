function [targets_est, stats] = analyze_results(g, targets, targets_est, name)

if ~exist('name','var')
    name = '';
end

P = perms(1:g.L);
hits = zeros(size(P));
goodness = zeros(size(P,1),1);
R = diag([g.hit_rate_threshold.t^-2 g.hit_rate_threshold.f^-2]);
%parfor p=1:size(P,1)
for p=1:size(P,1)
    t_error = targets.t - targets_est.t(P(p,:));
    overflow_errors = find(abs(t_error) > g.tau/2);
    t_error(overflow_errors) = t_error(overflow_errors) - sign(t_error(overflow_errors)) * g.tau;
%     assert(isempty(find(abs(t_error) > g.tau/2, 1)));
    
    f_error = targets.f - targets_est.f(P(p,:));
    overflow_errors = find(abs(f_error) > 1/g.tau/2);
    f_error(overflow_errors) = f_error(overflow_errors) - sign(f_error(overflow_errors)) * 1/g.tau;
%     assert(isempty(find(abs(f_error) > 1/g.tau/2, 1)));
    
    d = [t_error f_error].';
    D = diag(d'*R*d);
    hits(p,:) = D'<1;
    goodness(p) = sum(1./(D(D<1)));
end
best_perm_candidates = find(sum(hits,2) == max(sum(hits,2)));
[~, best_perm_index] = max(goodness(best_perm_candidates));
best_perm_index = best_perm_candidates(best_perm_index);
targets_est.best_perm = P(best_perm_index, :);
targets_est.best_perm_hits = hits(best_perm_index, :) > 0;
stats.num_hits = sum(targets_est.best_perm_hits);

%calculate errors
[~, best_perm_inv] = sort(targets_est.best_perm);

% time - fix errors for hits
t_est = targets_est.t(targets_est.best_perm);
t_error = targets.t - t_est;
overflow_hit_errors = find(abs(t_error) > g.tau/2 & targets_est.best_perm_hits');
t_est(overflow_hit_errors) = t_est(overflow_hit_errors) + sign(t_error(overflow_hit_errors)) * g.tau;
targets_est.t = t_est(best_perm_inv);
% time - fix errors for misses
t_error = targets.t - targets_est.t(targets_est.best_perm);
overflow_miss_errors = find(abs(t_error) > g.tau/2 & ~targets_est.best_perm_hits');
t_error(overflow_miss_errors) = t_error(overflow_miss_errors) - sign(t_error(overflow_miss_errors)) * g.tau;
% assert(~any(abs(t_error) > g.tau/2));

% frequency - fix errors for hits
f_est = targets_est.f(targets_est.best_perm);
f_error = targets.f - f_est;
overflow_hit_errors = find(abs(f_error) > 1/g.tau/2 & targets_est.best_perm_hits');
f_est(overflow_hit_errors) = f_est(overflow_hit_errors) + sign(f_error(overflow_hit_errors)) * 1/g.tau;
targets_est.f = f_est(best_perm_inv);
% frequency - fix errors for misses
f_error = targets.f - targets_est.f(targets_est.best_perm);
overflow_miss_errors = find(abs(f_error) > 1/g.tau/2 & ~targets_est.best_perm_hits');
f_error(overflow_miss_errors) = f_error(overflow_miss_errors) - sign(f_error(overflow_miss_errors)) * 1/g.tau;
% assert(~any(abs(f_error) > 1/g.tau/2));

misses = ~targets_est.best_perm_hits;
stats.num_misses.t_only = sum(abs(t_error(misses)) > g.hit_rate_threshold.t & ...
    abs(f_error(misses)) < g.hit_rate_threshold.f);
stats.num_misses.f_only = sum(abs(f_error(misses)) > g.hit_rate_threshold.f & ...
    abs(t_error(misses)) < g.hit_rate_threshold.t);
stats.num_misses.t_and_f = sum(abs(t_error(misses)) > g.hit_rate_threshold.t & ...
    abs(f_error(misses)) > g.hit_rate_threshold.f);

if stats.num_hits > 0
    stats.mse_t = mean(t_error(targets_est.best_perm_hits).^2);
    stats.mse_f = mean(f_error(targets_est.best_perm_hits).^2);
    stats.bias_t = mean(t_error(targets_est.best_perm_hits));
    stats.bias_f = mean(f_error(targets_est.best_perm_hits));
else
    stats.mse_t = nan;
    stats.mse_f = nan;
    stats.bias_t = nan;
    stats.bias_f = nan;
end

if ~isempty(name)
    fprintf('%s: hit_rate=%.2f RMS(t)=%.2f RMS(f)=%.2f bias(t)=%.2f bias(f)=%.2f num misses(t/f/both)=%d/%d/%d\n', ...
        name, mean(targets_est.best_perm_hits), sqrt(stats.mse_t)/g.Nyquist.dt, sqrt(stats.mse_f)/g.Nyquist.df, ...
        stats.bias_t/g.Nyquist.dt, stats.bias_f/g.Nyquist.df, stats.num_misses.t_only, ...
        stats.num_misses.f_only, stats.num_misses.t_and_f);
end