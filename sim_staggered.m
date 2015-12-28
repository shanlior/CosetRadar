function [successVec,resultHist,realHist,targets,targets_staggered] = sim_staggered(Ci,Q,L,P,snr_db,plot_fail_sim,numSims)
if ( nargin == 0) 
%    Ci=[0 3 5 7 11  17 19 23 ]; % channel coefficient
    Ci=[0];
    Q = 4; % How many ambiguities are resolved
    L = 5; % numTargets
%     P= 10; % numPulses Kron
    P = 100;
    plot_fail_sim = false;
    numSims = 30;
    snr_db = -20;
%     numSims = 1;
end
% Simulation config
rng('shuffle');
success = 0;
sumHits = 0;
realHist = zeros(numSims,L,2);
resultHist = zeros(numSims,L,2);
successVec = zeros(numSims,1);
rngVec=1:numSims;
failVec=[4,17,44,50,75,84,97,100];

for (i=1:numSims)
    disp(['Simulation number ',int2str(i)]);
    Results(i).a=rng(rngVec(i));
    Results(i).a=rng('shuffle');

    g_coset = global_settings(P,P,L, Ci,Q,snr_db);
    g_coset.numSims = numSims;

    % Randomizing targets
    tau = [1,1,1.05,1.05,1.1,1.1]*g_coset.tau;
%     tau = g_coset.tau;
    targets = randomize_targets(g_coset);
%     round(targets.t*g_coset.Fs)
%     round(g_coset.tau*g_coset.Fs)
    targets_staggered = MultiplePRF(tau, g_coset,targets);
    [isSuccess,realHist(i,:,:),resultHist(i,:,:),successVec(i)] = ...
        analyze_result(g_coset,targets,targets_staggered,i,plot_fail_sim,true);
     success = success + isSuccess;
    sumHits = sumHits + successVec(i);
    fprintf('Success rate: %.1f percent\n', 100*success/i);
    fprintf('Hit Rate: %.2f\n', sumHits/(i*L));
%     [targets_est, stats] = analyze_results_staggered(g_coset,...
%         targets, targets_staggered,'Staggered')
%     hit_rate = sum(abs(sort(targets.t) - sort(targets_staggered.t)) <= 3*g_coset.CS.delta_t);
%     sumHits = sumHits + stats.num_hits;
%     if stats.num_hits == L 
%         success = success + 1;
%     end
end

% sort(targets.t)
% sort(targets_staggered.t')
if success == numSims
    fprintf('Great Success!\n');
else
    disp([num2str(success), ' out of ' ,num2str(numSims)]);
    disp(['Hit Rate: ' num2str(sumHits/L)]);
end
end


% targets_Coset.f
% targets_Coset.t
% plot_results(g_coset, targets, targets_Coset);