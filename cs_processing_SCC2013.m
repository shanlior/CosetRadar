function [targets] = cs_processing_SCC2013(g, x)

% Xampling
C = zeros(length(g.CS.kappa), g.P);
X = fft(x);
kappa_indexes = mod(g.CS.kappa, size(X,1));
H_kappa = get_H_empiric(g);
for p=1:g.P
    C(:,p) = X(kappa_indexes+1,p);
    if g.CS.normalize_H_with_division
        C(:, p) = C(:, p) ./ H_kappa;
    else
        C(:, p) = C(:, p) .* H_kappa;
    end
end
clear X x;

if g.P > 1 && 0
    if 0
        [Omega, X_k] = OMP(g.CS.A2, C(:), g.L, 0, 1, g.CS.N_t);
    else
        [Omega, X_k] = OMP2(g, length(g.CS.kappa)*g.P, C(:), g.L, 1, g.CS.N_t);
    end
    targets.a = X_k;
    targets.t = mod(Omega-1,g.CS.N_t) / g.CS.N_t * g.tau;
    targets.f = floor((Omega-1)/g.CS.N_t) / g.CS.N_f / g.tau;
    return;
end

% CS recovery - first stage (delays)
L_delay = g.L * g.CS.time_over_recovery;
if 0
    [Omega, X_k, time_corrections] = OMP_Fourier(g.CS.A, g.CS.kappa, C, L_delay, g.CS.N_t, 3);
elseif 1
    [Omega, X_k, time_corrections] = OMP(g.CS.A, C, L_delay, 0, 3);
elseif 0
    addpath('Exponentials_toolbox');
    assert(all(kappa_indexes == sort(kappa_indexes)));
    [Omega, X_k] = MP_fb(conj(C), round(0.67*length(g.CS.kappa)), L_delay);
    Omega = round(Omega*g.CS.N_t) + 1;
    time_corrections = zeros(L_delay, 1);
elseif 0
    X = CoSaMP(g.CS.A, C, L_delay);
    Omega = find_nonzero_rows(X);
    X_k = X(Omega, :);
    time_corrections = zeros(L_delay, 1);
elseif 0
    if 1
        D_vec = L_delay:-1:1;
    else
        D_vec = ones(L_delay, 1);
        D_vec(1) = 3;
        D_vec(2) = 2;
    end
    X = MBMP_wrapper(g.CS.A, C, D_vec);
    Omega = find_nonzero_rows(X);
    X_k = X(Omega, :);
    time_corrections = zeros(L_delay, 1);
else
    X = RA_ORMP(g.CS.A, C, L_delay);
    Omega = find_nonzero_rows(X);
    X_k = X(Omega, :);
    time_corrections = zeros(L_delay, 1);
end

% if 0
%     [X] = gap_cs(g.CS.A, C, g.L);
% elseif 0
%     [wk,xk] = music_f(conj(C), g.L);
%     X = zeros(g.CS.N_t, g.P);
%     X(round(wk*length(X))) = xk;
% elseif 0
%     [x_L1] = l1l2(g.CS.A, C, L);
% end

assert(length(Omega) == L_delay);
assert(all(abs(time_corrections) <= 0.5));
CS_times = (Omega - 1 + time_corrections) * g.CS.delta_t;

% CS recovery - second stage (Doppler)
if g.P == 1
    L_freq = 1;
else
    L_freq = g.L;
end
targets_tmp.a = zeros(L_freq * L_delay, 1);
targets_tmp.t = zeros(L_freq * L_delay, 1);
targets_tmp.f = zeros(L_freq * L_delay, 1);
for k=1:L_delay
    if g.P == 1
        CS_freq = 0;
        x_q = X_k(k, 1);
    else
        a_q = X_k(k, :).';
        if 1
            [Omega, x_q, freq_corrections] = OMP(g.CS.B, a_q, L_freq, 1, 3);
        elseif 0
            [Omega, x_q, freq_corrections] = OMP_Fourier(g.CS.B, 0:g.P-1, a_q, L_freq);
        elseif 0
            MatrixPencil_KfirRonen(a_q,L_freq,T,nn,htilde);
        else
            A_q = abs(fft(a_q));
            x_q_fft = zeros(size(x_q));
            [best_freq_power, best_freq] = sort(A_q,'descend');
            x_q_fft(round(best_freq(1:L_freq)/length(a_q)*length(x_q))) = best_freq_power(1:L_freq);
            x_q = x_q_fft;
        end
        assert(length(Omega) == L_freq);
        assert(all(abs(freq_corrections) <= 0.5));
        CS_freq = (Omega - 1 + freq_corrections) * g.CS.delta_f;
    end
    
    targets_tmp.a(1+(k-1)*L_freq:k*L_freq) = x_q;
    targets_tmp.t(1+(k-1)*L_freq:k*L_freq) = CS_times(k);
    targets_tmp.f(1+(k-1)*L_freq:k*L_freq) = CS_freq;
end
if 1
    [~,strong_indexes] = sort(abs(targets_tmp.a), 'descend');
else
    tmp = sum(abs(X_k),2);
    [~,strong_indexes] = sort(tmp, 'descend');
end
targets.a = targets_tmp.a(strong_indexes(1:g.L));
targets.t = targets_tmp.t(strong_indexes(1:g.L));
targets.f = targets_tmp.f(strong_indexes(1:g.L));
