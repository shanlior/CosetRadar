function [targets] = randomize_targets(g)

% randomize amplitudes
targets.a = ones(g.L, 1);
if g.P > 1
    targets.a = targets.a .* exp(1j*2*pi*rand(g.L,1));
end
if ~g.fixed_target_amplitudes % add +/-10dB to amplitude, uniformly distributed in dB
    delta_a_max_dB = 10;
    delta_a = delta_a_max_dB * (2*rand(g.L, 1) - 1);
    targets.a = targets.a .* 10^(delta_a / 10);
end

% randomize time and frequency, make sure no overlap
overlapping_targets = 1;
min_separation = 3;    
while overlapping_targets
    % g.L = number of targets
    % t_pulse - the width of the pulse in the time domain. we need to
    % substract it because we don't want to create another bucket
    % we added the Q factor
    targets.t = rand(g.L, 1) * (g.Q * g.tau - g.t_pulse); %[sec]
    targets.t = round(targets.t * 1e8)/1e8;
    targets.f = rand(g.L, 1) * 1/g.tau; %[hz]
    overlapping_targets = 0;
    for l=1:g.L-1
        for ll=l+1:g.L
            % needs to check if seperation is limited because of tau
            % seperation
            if abs(targets.t(l) - targets.t(ll)) < min_separation*g.Nyquist.dt && ...
                    abs(targets.f(l) - targets.f(ll)) < min_separation*g.Nyquist.df
                overlapping_targets = 1;
                break;
            end
        end
    end
end


% error('start using target dipoles in order to really check resolution limits');