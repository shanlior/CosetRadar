function [Omega, X_k, Omega_corrections] = OMP(A, Y, K, do_interp, nn_erase, NN)
% implement an MMV version of OMP (S-OMP)
% inputs:
% A - sensing matrix, MxN
% Y - measurement vector, MxP
% k_max - max number of iterations to run
% do_interp - perform local interpolation
% nn_erase - number of nearest neighbours to erase around selected indexes
% outputs:
% Omega - active rows in X
% X_k - values in Omega
% Omega_corrections - off-grid corrections to Omega (if do_interp)

if ~exist('do_interp','var')
    do_interp = 0;
end
if ~exist('nn_erase','var')
    nn_erase = 0;
end

[M, N] = size(A);
assert(size(Y,1) == M);
assert(M > K);
Omega_corrections = zeros(K, 1);
AA = zeros(M, K);
R = Y;
Omega = []; % support
Q = [.5 -1 .5 ; -.5 0 .5 ; 0 1 0]; % for parabola interpolation
for k=1:K
    B = abs(A'*R);
    b = sum(B.^2, 2);
    
    b(Omega) = 0;
    for nn=1:nn_erase
        b(min(length(b),Omega+nn)) = 0;
        b(max(1,Omega-nn)) = 0;
        if exist('NN','var')
            b(min(length(b),Omega+NN*nn)) = 0;
            b(max(1,Omega-NN*nn)) = 0;
        end
    end
    
    [~, new_index] = max(b); % new_index = find(Thresholding(b, 1));
    Omega = [Omega ; new_index]; % update support
    assert(length(Omega) == k);
    
    if do_interp && new_index > 1 && new_index < N
        p = b(new_index+(-1:1));
        a = Q * p;
        assert(a(1) < 0);
        Omega_corrections(k) = -a(2) / 2 / a(1);
    end
    
    % update X estimate
    if Omega_corrections(k) == 0
        AA(:,k) = A(:, new_index);
    elseif Omega_corrections(k) > 0
        AA(:,k) = (1-Omega_corrections(k)) * A(:, new_index) + Omega_corrections(k) * A(:,new_index+1);
    else
        AA(:,k) = (1+Omega_corrections(k)) * A(:, new_index) - Omega_corrections(k) * A(:,new_index-1);
    end
    R = Y - AA(:,1:k) * pinv(AA(:,1:k)) * Y; % update residual
end
X_k = pinv(AA) * Y;
assert(all(abs(Omega_corrections) <= 0.5));
