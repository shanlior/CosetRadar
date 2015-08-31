function [targets] = coset_nyquist(g, x, SolveAlgorithm) % SolveAlgorithm=1=FISTA; 2=OMP

if nargin == 2
    SolveAlgorithm = 2; % Default value is OMP
end

N = g.CS.N_t;
P = g.P;
Ci=g.Ci;

k = g.CS.kappa;
m_p = g.m_p;

% Xampling
H_kappa = get_H_empiric(g);
C=zeros(length(g.CS.kappa), P*length(Ci));
for u=1:length(Ci)
    C_tmp = zeros(length(g.CS.kappa), P);
%     X = fft(x(:,1:size(C,2)));
    X = fft(x(:,((u-1)*size(C_tmp,2)+1):(u*size(C_tmp,2))));
    kappa_indexes = mod(g.CS.kappa, size(X,1));
    for p=1:size(X,2)
        C_tmp(:,p) = X(kappa_indexes+1,p);
        if g.CS.normalize_H_with_division
            C_tmp(:, p) = C_tmp(:, p) ./ H_kappa;
        else
            C_tmp(:, p) = C_tmp(:, p) .* H_kappa;
        end
        C(:,((u-1)*size(C_tmp,2)+1):(u*size(C_tmp,2)))=C_tmp;
    end
    
end
sizex=size(x)
sizeX=size(X)
clear X;

% Creates the Vandermonde matrices

% --------------------------------------------- Added
F=10; % how much ambiguities are resolved
M=P*F;
F_M_cP = zeros(M,P*length(Ci)); % for all channels

% TODO: check the 2*pi
F_M_P = (exp(-j*2*pi*((0:P*F-1)')*(0:P-1)/P));
for i = 1:length(Ci)
    tempMat = diag(exp(-j*2*pi*Ci(i)*(0:P-1)/P));
    F_M_cP(:,(1:P) + (i-1)*P) = F_M_P*tempMat;
end


F_K_N = exp(-j*2*pi*k*(0:N-1)/N);
N
sizek=size(k)
sizeC=size(C)
sizeC_tmp=size(C_tmp)
sizeF_K_N=size(F_K_N)
sizeF_M_cP=size(F_M_cP)


% Solve algorithm
if SolveAlgorithm==1
    disp('Running FISTA');
    map = FISTA(C,F_K_N,F_M_cP',N,M);
else
    disp('Running OMP');
    [map, Supp, x_t] = OMPmatrix(C,F_K_N,F_M_cP,20);
end

map_power = abs(map).^2;

% Peak detection
sample_SubNyquist_factor = 1;
[~,map_power_indexes] = sort(map_power(:),'descend');
targets.a = zeros(g.L, 1);
targets.t = zeros(g.L, 1);
targets.f = zeros(g.L, 1);
l = 0;
index_count = 1;
local_max = ordfilt2(map_power, 9, ones(3));
Q = [.5 -1 .5 ; -.5 0 .5 ; 0 1 0]; % for parabola interpolation
while l < g.L && index_count <= length(map_power_indexes)
    [t_index,f_index] = ind2sub(size(map_power), map_power_indexes(index_count));
    if local_max(t_index, f_index) == map_power(t_index, f_index)
        l = l + 1;
        targets.a(l) = map(t_index, f_index);
        time_bias_fix = 0.4;
        
        % time intepolation
        t_index_fix = 0;
        if t_index > 1 && t_index < size(map_power, 1)
            t_power = map_power(t_index+(-1:1), f_index);
            a = Q * t_power;
            assert(a(1) < 0);
            t_index_fix = -a(2) / 2 / a(1);
        end
        targets.t(l) = ((t_index - 1 + t_index_fix) * sample_SubNyquist_factor + time_bias_fix) * g.tau/N;
        assert(targets.t(l)>=0 && targets.t(l)<=g.tau);
        
        % frequency interpolation
        f_index_fix = 0;
        if f_index > 1 && f_index < size(map_power, 2)
            f_power = map_power(t_index, f_index+(-1:1));
            a = Q * f_power.';
            assert(a(1) < 0);
            f_index_fix = -a(2) / 2 / a(1);
        end
        targets.f(l) = (f_index - 1 + f_index_fix) /(M*g.tau);
        assert(targets.f(l)>=0 && targets.f(l)<=1/g.tau);
        
        if g.mf.debug_plot
            h = plot(targets.t(l), targets.f(l), 'ro');
            set(h,'MarkerSize',10);
        end
    end
    index_count = index_count + 1;
end
