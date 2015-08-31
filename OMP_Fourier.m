function [Omega, X_k, Omega_corrections] = OMP_Fourier(A, kappa, Y, K, N_t, nn_erase)
% implement an MMV version of OMP (S-OMP)
% inputs:
% A - sensing matrix, MxN
% Kappa - vector of frequency bins
% Y - measurement vector, MxP
% K - sparsity level
% nn_erase - number of nearest neighbours to erase around selected indexes
% outputs:
% Omega - active rows in X
% X_k - values in Omega
% Omega_corrections - off-grid corrections to Omega

if ~exist('nn_erase','var')
    nn_erase = 0;
end

[M, N] = size(A);
assert(size(Y,1) == M);
assert(M > K);
assert(size(kappa,1) == M);
assert(isvector(kappa));
Omega_corrections = zeros(K, 1);
AA = zeros(M, K);
q = 2;
R = Y;
Omega = []; % support
Q = [.5 -1 .5 ; -.5 0 .5 ; 0 1 0]; % for parabola interpolation
for k=1:K
    B = abs(A'*R);
    b = sum(B.^q, 2).^(1/q);

    b(Omega) = 0;
    for nn=1:nn_erase
        b(min(length(b),Omega+nn)) = 0;
        b(max(1,Omega-nn)) = 0;
    end
    
    [~, new_index] = max(b); % new_index = find(Thresholding(b, 1));
    Omega = [Omega ; new_index]; % update support
    assert(length(Omega) == k);
    
    if new_index > 1 && new_index < N
        p = b(new_index+(-1:1)).^2;
        a = Q * p;
        assert(a(1) < 0);
        Omega_corrections(k) = -a(2) / 2 / a(1);
    end
    
    % update X estimate
    n = new_index - 1 + Omega_corrections(k);
    AA(:,k) = exp(-1j*2*pi*n*kappa/N_t);
    
    X_temp = pinv(AA(:,1:k)) * Y;
    R = Y - AA(:,1:k) * X_temp; % update residual
end
X_k = pinv(AA) * Y;
assert(all(abs(Omega_corrections) <= 0.5 + 1e-5));
