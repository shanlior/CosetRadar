function [x] = generate_analog_input_signal(g_coset, targets)

% if exist('g_coset.m_p','var')==1
%     P = size(g_coset.m_p,1);
% else
%     P = g_coset.P;
% end

P = g_coset.P;
Q = g_coset.Q;

% create "analog" signal
% buckets before chopping
num_buckets = P + Q - 1;
x = zeros(length(g_coset.Ci), round(g_coset.tau * g_coset.Fs), num_buckets);

% x_iter = zeros(round(g_coset.tau * g_coset.Fs), num_buckets); % for more than one channel we need to add a dimension
x_tmp = zeros(round(g_coset.tau * g_coset.Fs)*num_buckets,1); % a matrix for each channel (needed for another dimension support)
% Following from the equation: g_coset.tau * g_coset.Fs = Nsamples
%n_start = round(targets.t * g_coset.Fs) + 1;
n_start = round(targets.t * g_coset.Fs) + 1;
n_stop = n_start + g_coset.h_length - 1;
phi1 = -2j*pi*targets.f*(0:P-1)*g_coset.tau; % phase shift for calculations
phi2 = -2j*pi*targets.f.*floor(targets.t/g_coset.tau)*g_coset.tau; % phase shift for calculations

%===============================
C_delay_element=exp(-1j*2*pi*g_coset.Ci'*(0:P-1)/P) ;    % Phase shift by modulation
%===============================
for c=1:length(g_coset.Ci)
    for l=1:g_coset.L      
        for p=1:P
            x_tmp((n_start(l):n_stop(l)) + (p - 1) * round(g_coset.tau * g_coset.Fs),1) = ...
             x_tmp((n_start(l):n_stop(l)) + (p - 1) * round(g_coset.tau * g_coset.Fs),1) + ...
              targets.a(l) * g_coset.h * exp(phi1(l,p)+phi2(l))*C_delay_element(c,p);
        end
    end
    % Make a matrix out of the linearized signal
    
    x_iter = reshape(x_tmp,[round(g_coset.tau * g_coset.Fs), num_buckets]);  
    x(c,:,:) = x_iter;
end
    % chop matrix
    x = x(:,:,Q:P); % Dimensions: 1=Channel, 2=Time, 3=Bucket
    
    
    % add noise
    if g_coset.snr < inf % add noise
        Ps = get_h_power(g_coset);
        sigma_n = sqrt(Ps / g_coset.snr);
%         n = sigma_n * randn(size(x));
        n = sigma_n * crandn(size(x)) / sqrt(2);
        if 0
            figure;
            plot(real(n(1,:,1)),'g.:');
            hold on;
            plot(real(x(1,:,1)),'b.:');
        end
        x = x + n;
        
    end
end

function [Ps] = get_h_power(g)
Ps = g.h'*g.h / g.h_length; % get average power
end

