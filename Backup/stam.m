function [x] = stam(g, targets)

% if exist('g.m_p','var')==1
%     P = size(g.m_p,1);
% else
%     P = g.P;
% end

P = g.P;

% create "analog" signal
x = zeros(round(g.tau * g.Fs), P);
n_start = round(targets.t * g.Fs) + 1;
n_stop = n_start + g.h_length - 1;
phi = 2*pi*targets.f*(0:P-1)*g.tau;
for l=1:g.L
    for p=1:P
        x(n_start(l):n_stop(l), p) = x(n_start(l):n_stop(l), p) + ...
            targets.a(l) * g.h * exp(1j*phi(l,p));
    end
end

if g.snr < inf % add noise
    Ps = get_h_power(g);
    sigma_n = sqrt(Ps / g.snr);
    if P == 1
        n = sigma_n * randn(size(x));
    else
        n = sigma_n * crandn(size(x)) / sqrt(2); % sqrt(2) because noise is complex
    end
    if 0
        figure;
        plot(real(n(:,1)),'g.:');
        hold on;
        plot(real(x(:,1)),'b.:');
    end
    x = x + n;
end

if P == 1
    assert(all(isreal(x)));
end

function [Ps] = get_h_power(g)
Ps = g.h'*g.h / g.h_length; % get average power
