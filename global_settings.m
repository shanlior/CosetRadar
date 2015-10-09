function [g] = global_settings(nu_pulses, P, L, Ci,Q)

fftw('planner', 'hybrid');
g.L = L; % number of targets
g.P = P; % number of pulses
g.tau = 10e-6; % PRI [sec]
% g.tau = 5e-6; % PRI [sec]
g.h_BW = 1e6*11/10;
% g.h_BW = 100e6*2;
multFactor = 50;
g.Fs = multFactor*g.h_BW;
if 0
    addpath('OMP_Receiver_Test\SignalGeneration');
    addpath('OMP_Receiver_Test\Filters');
    [h, H_spectra] = genPulse(g.h_BW, g.Fs, g.tau);
    save('pulse','h','H_spectra');
else
    if 0
        load pulse;
    else
        h = zeros(15,1);
%         h(9) = -0.9;
%         h(7) = -0.9;
%         h(6) = 0.3;
%         h(10) = 0.3;
%         h(11) = -0.07;
%         h(5) = -0.07;
        h(8) = 3;
        H_spectra = fft(h,length(h)*multFactor);
    end
    
end
g.h = h;
g.h_length = length(h);
g.H_spectra = H_spectra;
g.t_pulse = g.h_length / g.Fs;
g.h_type = 2;
g.Nyquist.dt = 1/(2*g.h_BW*1); % Nyquist time bin [sec]
g.Nyquist.df = 1/(g.P*g.tau); % Nyquist frequency bin [hz]
g.hit_rate_threshold.t = 3*g.Nyquist.dt; % hit rate ellipse half time axis
g.hit_rate_threshold.f = 3*g.Nyquist.df; % hit rate ellipse half frequency axis
% g.snr = 10^(snr_db/10); %[1]
% g.fixed_target_amplitudes = fixed_target_amplitudes;
g.mf.use_windows = 0;
g.mf.debug_plot = 0;

% CS constants
g.CS.delta_t = g.Nyquist.dt/2;
g.CS.delta_f = g.Nyquist.df/2;
g.CS.N_t = round(g.tau / g.CS.delta_t);
g.CS.N_f = round(1/g.tau / g.CS.delta_f);
assert(abs(g.CS.N_t - g.tau / g.CS.delta_t) < 1e-11);
assert(abs(g.CS.N_f - 1/g.tau / g.CS.delta_f) < 1e-11);
g.invalid_indexes_end = floor(g.t_pulse/g.CS.delta_t)-1; % -1 to make even
g.CS.time_over_recovery = 2; % time over-recovery factor
g.CS.normalize_H_with_division = 1;

% assert(sample_SubNyquist_factor >= 1);
% assert(pulse_SubNyquist_factor >= 1);
k_max = round(0.95*g.h_BW*g.tau);
num_fourier_coeffs = round(2*g.h_BW*g.tau );
num_fourier_coeffs = min(num_fourier_coeffs, 2*k_max);

    g.CS.kappa = ceil(-num_fourier_coeffs/2):ceil(num_fourier_coeffs/2)-1;
    g.CS.kappa = 0:g.CS.N_t-1;
g.CS.kappa = g.CS.kappa(:);
% assert(length(g.CS.kappa) == num_fourier_coeffs);
% assert(all(abs(g.CS.kappa) <= k_max));

g.CS.A = single(get_V(g));
g.CS.B = get_B(g);
g.CS.A = normalize_columns(g.CS.A);
g.CS.B = normalize_columns(g.CS.B);

% g.CS.A2 = zeros(length(g.CS.kappa)*g.P, g.CS.N_t*g.CS.N_f, 'single');
% for p=0:g.P-1
%     g.CS.A2(1+p*length(g.CS.kappa):(p+1)*length(g.CS.kappa),:) = ...
%         kron(exp(1j*2*pi*p*(0:g.CS.N_f-1)/g.CS.N_f), g.CS.A);
% end
% assert(all(abs(g.CS.A2(:)) > 0));
if 0
    fprintf('duty cycle = %%%.2f\n', 100*g.t_pulse/g.tau);
    fprintf('Hit rate ellipse: +/-(%.2f,%.2f) Nyquist (time,freq) bins\n', ...
        g.hit_rate_threshold.t / g.Nyquist.dt, g.hit_rate_threshold.f / g.Nyquist.df);
    fprintf('CS grid bin over Nyquist bin ratio (time/freq) = (%.2f,%.2f)\n', ...
        g.CS.delta_t/g.Nyquist.dt, g.CS.delta_f/g.Nyquist.df);
    fprintf('Classic/CS processing using %d/%d samples per pulse, %d/%d pulses\n', ...
        round(g.tau*2*g.h_BW), length(g.CS.kappa), g.P, round(g.P / g.pulse_SubNyquist_factor));
    fprintf('Total CS undersampling ratio = %.1f\n', g.pulse_SubNyquist_factor * g.sample_SubNyquist_factor);
end

%% Gal Lior setting
% Adds m_p
if nargin == 0
    nu_pulses = 50;
end

g.m_p = sort(randsample(g.P,nu_pulses))-1;
% g.m_p = 1;
g.Ci = Ci;
g.Q = Q;