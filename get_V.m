function [V] = get_V(g)
assert(abs(g.CS.delta_t*g.CS.N_t-g.tau) < 1e-15);
V = exp(-1j*2*pi*g.CS.kappa*(0:g.CS.N_t-g.invalid_indexes_end-1)/g.CS.N_t);