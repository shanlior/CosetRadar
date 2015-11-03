clear all; clc;

N = 4;
M = 4;
m_p = [0,3]';
k = [3]';

A = exp(-j*2*pi*k*(0:N-1)/N);
B = (exp(-j*2*pi*m_p*(0:M-1)/M))';
X = zeros(4,4)*0.01;
X(2,3) = 10;
X(4,4) = 20;

Y = A*X*B;

[RES, SUPP, XT] = OMPmatrix(Y,A,B,2)

Y_vec = reshape(Y,length(k)*length(m_p),1);
X_vec = reshape(X,M*N,1);
B_kron_A = kron(transpose(B), A);

[Xr, Supp, ResNorm, NormResvsSol] = OMP_MMV(Y_vec, B_kron_A, 4, 0, 0);