function [g] = global_settings_nu(nu_pulses, P, L, snr_db, sample_SubNyquist_factor, pulse_SubNyquist_factor, ...
    random_fourier_coeffs, fixed_target_amplitudes,Ci)

% Builds default global settings
if nargin == 0 || nargin == 1
    g = global_settings();
elseif nargin == 2
    g = global_settings(P);
else
    g = global_settings(P, L, snr_db, sample_SubNyquist_factor, pulse_SubNyquist_factor, ...
    random_fourier_coeffs, fixed_target_amplitudes);
end

% Adds m_p
if nargin == 0
    nu_pulses = 50;
end

g.m_p = sort(randsample(g.P,nu_pulses))-1;
% g.m_p = 1;
g.Ci=Ci;
