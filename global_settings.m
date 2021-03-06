function [g] = global_settings(nu_pulses, P, L, Ci,Q, snr_db,...
    is_full_sample,reduce_method,sample_SubNyquist_factor,less_p,same_reduce_pulses_B)

fftw('planner', 'hybrid');
g.L = L; % number of targets
g.P = P; % number of pulses
g.tau = 10e-6; % PRI [sec]
% g.h_BW = 1e6*11/10; % kron
g.h_BW = 100e6; % normal
multFactor = 50;
g.Fs = multFactor*g.h_BW;
g.Ci = Ci;
g.Q = Q;
g.m_p = sort(randsample(g.P,nu_pulses))';
g.nu_pulses = nu_pulses;
g.same_reduce_pulses_B = 1;


random_fourier_coeffs = reduce_method;
if same_reduce_pulses_B 
    g.less_p = g.m_p;
else
    g.less_p =  sort(randsample(g.P,less_p))'; 
end
if 0
    addpath('OMP_Receiver_Test\SignalGeneration');
    addpath('OMP_Receiver_Test\Filters');
    [h, H_spectra] = genPulse(g.h_BW, g.Fs, g.tau);
    save('pulse','h','H_spectra');
else
    if 1
        assert(g.h_BW == 100e6);
        assert(g.Fs == 5e9);
        assert(g.tau == 10e-6);
        load pulse;
    else
        h = zeros(11,1);
% %         h(5) = -0.9;
%         h(7) = -0.9;
        h(10) = 0.05;
        h(8) = 0.5;
        h(4) = 0.5;
        h(2) = 0.05;
%         h(3) = -0.07;
%         h(9) = -0.07;
        h(6) = 3;
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
g.snr_db = snr_db;
g.snr = 10^(snr_db/10); %[1]
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


if is_full_sample == 0

    g.sample_SubNyquist_factor = sample_SubNyquist_factor;
    k_max = round(g.h_BW*g.tau);
    num_fourier_coeffs = round(2*g.h_BW*g.tau / g.sample_SubNyquist_factor);
    num_fourier_coeffs = min(num_fourier_coeffs, 2*k_max);
    if random_fourier_coeffs == 0 % LPF
        if g.P==1 && num_fourier_coeffs <= k_max
            disp('todo - real coeffs');
            g.CS.kappa = 0:num_fourier_coeffs-1;
        else
            g.CS.kappa = ceil(-num_fourier_coeffs/2):ceil(num_fourier_coeffs/2)-1;
        end
    elseif random_fourier_coeffs == 1 % true rand
        if g.P==1 && num_fourier_coeffs <= k_max+1
            g.CS.kappa = randsample(k_max+1, num_fourier_coeffs)-1;
        else
            g.CS.kappa = randsample(2*k_max+1, num_fourier_coeffs)-(k_max+1);
        end
    elseif random_fourier_coeffs == 2 % blocks of K
        K = 25;
        g.CS.kappa = randsample(2*k_max/K, num_fourier_coeffs/K)-1;
        g.CS.kappa = bsxfun(@plus, K*g.CS.kappa, 0:K-1);
        g.CS.kappa = g.CS.kappa(:) - k_max;
        assert(length(g.CS.kappa) == length(unique(g.CS.kappa)));
    elseif random_fourier_coeffs == 3 % cosets 1 - same k in all cosets
        if g.P==1
            k=floor(k_max/num_fourier_coeffs);
            m=randi(k)-1;
            g.CS.kappa=(m:k:k_max-1);
        else
            k=floor(2*k_max/num_fourier_coeffs);
            m=randi(k)-1;
            g.CS.kappa=(-k_max+m:k:k_max-1);
        end
    elseif random_fourier_coeffs == 4 % cosets 2 - different k in all cosets
        if g.P==1
            k=floor(k_max/num_fourier_coeffs);
            g.CS.kappa=randi(k, 1, num_fourier_coeffs)-1;
            g.CS.kappa=g.CS.kappa+(0:num_fourier_coeffs-1)*k;
        else
            k=floor(2*k_max/num_fourier_coeffs);
            g.CS.kappa=randi(k, 1, num_fourier_coeffs)-1;
            g.CS.kappa=g.CS.kappa+(0:num_fourier_coeffs-1)*k-k_max;
        end
    elseif random_fourier_coeffs == 5
        k_range=num_fourier_coeffs+200;
        g.CS.kappa=randsample(k_range,num_fourier_coeffs)-k_range/2;
    else
        error('bad option');
    end
    g.CS.kappa = g.CS.kappa(:);
else
    k_max = round(g.h_BW*g.tau);
    num_fourier_coeffs = round(4*g.h_BW*g.tau );
    g.CS.kappa = ceil(-num_fourier_coeffs/2):ceil(num_fourier_coeffs/2)-1;
    g.CS.kappa = g.CS.kappa(:);
end


if 0
    fprintf('duty cycle = %%%.2f\n', 100*g.t_pulse/g.tau);
    fprintf('Hit rate ellipse: +/-(%.2f,%.2f) Nyquist (time,freq) bins\n', ...
        g.hit_rate_threshold.t / g.Nyquist.dt, g.hit_rate_threshold.f / g.Nyquist.df);
    fprintf('CS grid bin over Nyquist bin ratio (time/freq) = (%.2f,%.2f)\n', ...
        g.CS.delta_t/g.Nyquist.dt, g.CS.delta_f/g.Nyquist.df);
    fprintf('Classic/CS processing using %d/%d samples per pulse, %d/%d pulses\n', ...
        round(g.tau*2*g.h_BW), length(g.CS.kappa), g.P, g.nu_pulses);
    fprintf('Total CS undersampling ratio = %.1f\n', (g.P/g.nu_pulses) * g.sample_SubNyquist_factor);
end


