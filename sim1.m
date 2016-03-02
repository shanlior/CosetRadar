function [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,snr_db,plot_fail_sim,numSims,is_full_sample,reduce_method,...
    sample_SubNyquist_factor,nu_pulses,less_p,same_reduce_pulses_B)
if ( nargin == 0) 
%    Ci=[0 3 5 7 11  17 19 23 ]; % channel coefficient
    Ci=[0 3 5 7 11 13 17];
    Q = 4; % How many ambiguities are resolved
    L = 5; % numTargets
%     P= 10; % numPulses Kron
    P = 100;
    plot_fail_sim = false;
    numSims = 5;
    snr_db = -30;
    snr_db = inf;
    is_full_sample = 1;
    reduce_method = 1;
    sample_SubNyquist_factor = 1;
    nu_pulses = P;
%     numSims = 1;
    less_p = P;
    same_reduce_pulses_B = 1;
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
    Results(i).a=rng(rngVec(i));
    Results(i).a=rng('shuffle');
    
    g_coset = global_settings(nu_pulses,P,L, Ci,Q,snr_db,is_full_sample,...
        reduce_method,sample_SubNyquist_factor,less_p,same_reduce_pulses_B);
    g_coset.numSims = numSims;

    % Randomizing targets
    targets = randomize_targets(g_coset);

    % Generating radar received analog signal (slow-fast matrix)
    x = generate_analog_input_signal(g_coset, targets);
    % x = stam(g, targets);
%  disp {coset_nyquist start}
    % Processing
    targets_Coset = coset_nyquist(g_coset,x,targets);
%  disp {coset_nyquist end}
%     if g.mf.debug_plot
%         for l=1:g.L
%             plot(targets.t(l), targets.f(l), 'bo');
%         end
%     end

    [isSuccess,realHist(i,:,:),resultHist(i,:,:),successVec(i)] = analyze_result(g_coset,targets,targets_Coset,i,plot_fail_sim);
    success = success + isSuccess;
    sumHits = sumHits + successVec(i);
    fprintf('Success rate: %.1f percent\n', 100*success/i);
    fprintf('Hit Rate: %.2f\n', sumHits/(i*L));
   
end
if success == numSims
    fprintf('Great Success!\n');
end

% targets_Coset.f
% targets_Coset.t
% plot_results(g_coset, targets, targets_Coset);