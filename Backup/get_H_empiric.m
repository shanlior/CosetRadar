function [H] = get_H_empiric(g)

if g.h_type == 0 || g.h_type == 1
    k_Nyquist = round(g.Nyquist.f * g.tau / 2);
    assert(all(abs(g.CS.kappa) <= k_Nyquist));
    
    F = g.tau / g.t_pulse;
    h = get_h(g, 0, 0:g.Nyquist.dt:g.t_pulse-g.Nyquist.dt);
    H_all = g.t_pulse / length(h) * fftshift(fft(h,1+F*length(h)));
    H = H_all(g.CS.kappa + k_Nyquist + 1);
elseif g.h_type == 2
    kappa_indexes = mod(g.CS.kappa, length(g.H_spectra));
    H = g.H_spectra(kappa_indexes + 1);
else
    error('no such type');
end