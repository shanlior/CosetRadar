function [targets] = randomize_targets(g)

% randomize amplitudes
targets.a = ones(g.L, 1);
if g.P > 1
    targets.a = targets.a.* exp(1j*2*pi*rand(g.L,1));
end

% randomize time and frequency, make sure no overlap
overlapping_targets = 1;
min_separation = 3;    
while overlapping_targets
    % t_pulse - the width of the pulse in the time domain. we need to
    % substract it because we don't want to create another bucket
    % we added the Q factor
    targets.t = rand(g.L, 1) * (g.Q * g.tau - g.t_pulse); %[sec]
    targets.t = round(targets.t * g.Fs) / g.Fs;
    targets.t = floor(targets.t / g.CS.delta_t) * g.CS.delta_t;
    targets.f = rand(g.L, 1) * 1/g.tau; %[hz]
    targets.f = floor(targets.f *  g.P * g.tau) / (g.P * g.tau); 
    overlapping_targets = 0;
    for l=1:g.L-1
        for ll=l+1:g.L
            if abs(targets.t(l) - targets.t(ll)) < min_separation*g.CS.delta_t
                overlapping_targets = 1;
                break;
            end
        end
    end
end
