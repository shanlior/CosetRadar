function [ targets_staggered ] = MultiplePRF( tau, g,targets )
% Solving radar signal using the multiplePRF scheme
    debug_print = false;
    P = g.P;
    L = g.L;
    Ci = g.Ci;
    Q = g.Q;
    snr_db = g.snr_db;
    LCM = 1;
    for c=1:length(tau)
        LCM = lcm(LCM,round(tau(c)*g.Fs));
    end
    for c=1:length(tau)
        if debug_print
            disp(['Processing channel ',int2str(c)]);
        end
        g_new = global_settings_tau(P,P,L, Ci,Q,snr_db,tau(c));
        x{c} = generate_analog_input_signal_staggered(tau(c),g_new, targets);
        g_new = global_settings_tau(P-Q+1,P-Q+1,L, Ci,Q,snr_db,tau(c));
        targets_tmp = classic_processing(g_new, x{c});
        targets_ch(c) = targets_tmp;
        
%         figure;
%         plot(abs(x{c}))
%         [targets_ch(c).t,ind]=sort(targets_ch(c).t);
%         targets_ch(c).f=targets_ch(c).f(ind);
%         targets_ch(c).a=targets_ch(c).a(ind);
    end
    
    if debug_print
        disp(['Disecting results']);
    end

    t = zeros(length(tau),LCM);
    for l=1:L
        for c=1:length(tau)
            Qsize = LCM / round(tau(c)*g.Fs);
            for q=0:Qsize-1
                indices = (round(targets_ch(c).t(l) * g.Fs) + 1 + q *...
                        round(tau(c)*g.Fs))+ ...
                        ((ceil(-g.h_length/2)+1):1:ceil(g.h_length/2));
                indices = indices(find(indices > 0));
                firstIdx = length(g.h) - length(indices) + 1;
                t(c,indices) = g.h(firstIdx:end);

            end
        end
    end
    targets_tmp = ones(1,LCM);
    for c=1:size(t,1)
        targets_tmp = targets_tmp .* t (c,:);
    end
    numPeaks = 0;
    peakDist = g.h_length;
    while (numPeaks ~= 5) % making sure 5 targets are found
        [~,locs] = findpeaks(targets_tmp,'SORTSTR','descend','NPEAKS',L ,...
            'MINPEAKDISTANCE',peakDist);
        numPeaks = length(locs);
        peakDist = floor(peakDist / 2);
        if ~peakDist
            break;
        end
    end
    targets_staggered.t = locs / g.Fs;
    
    

    
    % finding t - q=0
    targets_q0 = rem(targets_staggered.t,tau(1));
    for l=1:length(targets_q0)
        [~,targetIdx] = min(abs(targets_ch(1).t - targets_q0(l)));
        targets_staggered.f(l) = targets_ch(1).f(targetIdx);
        targets_staggered.a(l) = targets_ch(1).a(targetIdx);
    end
    

    targets_staggered.t = targets_staggered.t.';
    targets_staggered.f = targets_staggered.f.';
    targets_staggered.a = targets_staggered.a.';

    
end

