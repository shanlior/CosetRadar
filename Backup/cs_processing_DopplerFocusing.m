function [targets] = cs_processing_DopplerFocusing(g, x)

% Xampling
C = zeros(length(g.CS.kappa), round(g.P / g.pulse_SubNyquist_factor));
X = fft(x(:,1:size(C,2)));
kappa_indexes = mod(g.CS.kappa, size(X,1));
H_kappa = get_H_empiric(g);
for p=1:size(X,2)
    C(:,p) = X(kappa_indexes+1,p);
    if g.CS.normalize_H_with_division
        C(:, p) = C(:, p) ./ H_kappa;
    else
        C(:, p) = C(:, p) .* H_kappa;
    end
end
clear X x;

if 0
    D=C;
    load C;
    figure;
    hist(20*log10(abs([C(:) D(:)])),1000);
    xlabel('|c_p[k]|^2 [dB]');
    ylabel('# occurrences');
    legend('target','clutter');
end

if 0
    w = taylorwin(size(C,2),4,-50).';
    % w = chebwin(size(C,2), 50).';
    C = C .* repmat(w,size(C,1),1);
end

targets.a = zeros(g.L,1);
targets.t = zeros(g.L,1);
targets.f = zeros(g.L,1);
for l=1:g.L
    Psi = fft(C, g.CS.N_f, 2); % Doppler focusing
    
    if 0
        K = 5;
        Psi(:,1:K) = 0;
        Psi(:,end-(K-1):end) = 0;
    end
    
    for f=1:g.CS.N_f
        [Omega, X_k, time_corrections] = OMP_Fourier(g.CS.A, g.CS.kappa, Psi(:,f), g.L, g.CS.N_t, 3);
        assert(length(Omega) == g.L);
        assert(all(abs(time_corrections) <= 0.5 + 1e-4));
        
        if max(abs(X_k)) > abs(targets.a(l))
            [~, best_index] = max(abs(X_k));
            targets.a(l) = X_k(best_index);
            targets.t(l) = (Omega(best_index) - 1 + time_corrections(best_index)) * g.CS.delta_t;
            targets.f(l) = (f-1) / g.CS.N_f / g.tau;
        end
    end
    
    % subtract estimated target response
    targets.a(l) = targets.a(l) / g.P;
    
    C_temp = targets.a(l) * exp(-1j*2*pi*g.CS.kappa*targets.t(l)/g.tau) * ...
        exp(1j*2*pi*targets.f(l)*(0:g.P-1)*g.tau);
    C = C - C_temp;
end