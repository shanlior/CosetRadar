function [B] = get_B(g)
assert(abs(g.CS.delta_f*g.CS.N_f-1/g.tau) < 1e-10);
B = exp(1j*2*pi*(0:g.P-1).'*(0:g.CS.N_f-1)/g.CS.N_f);