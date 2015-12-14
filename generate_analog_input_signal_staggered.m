function [x] = generate_analog_input_signal_staggered(tau,g_coset,targets);

    P = g_coset.P;
    Q = g_coset.Q;

    % create "analog" signal
    % buckets before chopping
    num_buckets = P + Q - 1;
    x = zeros(round(tau * g_coset.Fs) * num_buckets,1);

    % Following from the equation: g_coset.tau * g_coset.Fs = Nsamples



    % % create transmitted signal
    % 
    % [~,pulseLocation] = max(g_coset.h);
    % for c=1:numChannels   
    %     for p=1:P         
    %         x(c,pulseLocation + (p - 1) * round(tau * g_coset.Fs)) = ...
    %          x(c,pulseLocation + (p - 1) * round(tau * g_coset.Fs)) + ...
    %           1;
    %     end
    % end


    n_start = round(targets.t * g_coset.Fs) + 1;
    n_stop = n_start + g_coset.h_length - 1;

    % create recieved signal

    phi1 = 2j*pi*targets.f*(0:P-1)*tau; % phase shift for calculations
%     phi2 = -2j*pi*targets.f.*floor(targets.t/tau)*tau; % phase shift for
%     calculations CHECK
    for l=1:g_coset.L      
        for p=1:P     
            x((n_start(l):n_stop(l)) + (p - 1) * round(tau * g_coset.Fs)) = ...
             x((n_start(l):n_stop(l)) + (p - 1) * round(tau * g_coset.Fs)) + ...
             targets.a(l) * g_coset.h * exp(phi1(l,p));
%                targets.a(l) * g_coset.h * exp(phi1(l,p)+phi2(l)); CHECK

        end
    end
    
    x = reshape(x,[round(tau * g_coset.Fs), num_buckets]);  %  1=Time, 2=Bucket

    % chop matrix
    x = x(:,Q:P);

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
    
%     
%     % decimate to Nyquist rate
%     Li = round(g_coset.Fs/(2*g_coset.h_BW));
%     h_Nyquist = decimate(g_coset.h,Li);
% 
%     x_Nyquist = zeros(size(x,1),ceil(size(x,2)/Li));
%     for c=1:size(x,1)
%         x_Nyquist(c,:) = decimate(x(c,:),Li);
%         mapTmp(c,:) = xcorr(x_Nyquist(c,:), h_Nyquist);
%     end
%     map = mapTmp(:,length(h_Nyquist):end-length(h_Nyquist));
%     map_power = abs(map).^2;
%     
% %   find peaks
%     for c = 1:size(x,1)
%         [~,locs(c,:)] = findpeaks(map_power(c,:),'MINPEAKDISTANCE',floor(length(h_Nyquist)/2),...
%            'SORTSTR','descend','NPEAKS',P);
%     end
% 
% %     chop matrix
%     map_power = map_power(:,round(LCM/Li):end); % Dimensions: 1=Channel, 2=Time
%     
%     % chop locations
%     locs = locs - floor(length(h_Nyquist)/2) - LCM/Li;
% %     for c=1:size(locs,1)
% %         for i=1:size(locs,2)
% %             locs(c,i) = 1 - LCM + (locs(c,i))*Li
% %         end
% %     end
%     locs = sort(locs,2);
%     for c = 1:size(x,1)
%        idx(c) = min(find(locs(c,:)>0))
%     end    
%     
%      figure
%     plot(map_power(1,:))
%     
%     
% 
%    
% %     [~,map_power_indexes] = sort(map_power,2,'descend');
% %     local_max = zeros(size(map));
% %     % ordfilt 1 dimension
% %     for c=1:size(map,1)
% %         for i=2:size(map,2)-1
% %             local_max(c,i) = max(map_power(c,i-1:i+1));
% %         end
% %         local_max(c,1) = max(map_power(c,1:2));
% %         local_max(c,end) = max(map_power(c,(end-1):end));
% %     end
% %     
% %     l=0;
% %     index_count = 0;
% %     while l < g_coset.L && index_count <= size(map_power_indices,2)
% %         if local_max(
% %         
% %     end
%     
% 
%     % find peaks without noise
% %      for c = 1:size(x,1)
% %        [~,locs_noiseless(c,:)] = findpeaks(abs(x(c,:)),'MINPEAKDISTANCE',ceil(length(g_coset.h)/2),...
% %         'SORTSTR','descend','NPEAKS',P);
% %      end
% %      locs_noiseless = locs_noiseless - ceil(LCM/Li) - floor(g_coset.h_length/2);
% %      locs_noiseless = sort(locs_noiseless,2)
% 
% 
%     
%     

        
end

function [Ps] = get_h_power(g)
Ps = g.h'*g.h / g.h_length; % get average power
end

