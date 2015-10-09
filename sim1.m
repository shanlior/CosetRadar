function [successVec,resultHist,realHist,Results,targets,targets_Coset] = sim1(Ci,Q,L,P)
close all;
if ( nargin == 0) 
    Ci=[0 3 7 8 11 13 15 17 19]; % channel coefficient
    Q = 4; % How many ambiguities are resolved
    L = 4; % numTargets
    P=10; % numPulses
end
% Simulation config
rng('shuffle');
numSims = 100;
numSims = 20;
success = 0;
realHist = zeros(numSims,L,2);
resultHist = zeros(numSims,L,2);
successVec = zeros(numSims,1);
rngVec=1:numSims;
failVec=[4,17,44,50,75,84,97,100];

for (i=1:numSims)

%     Results(i).a=rng(rngVec(i));
    Results(i).a=rng('shuffle');

    g_coset = global_settings(P,P,L, Ci,Q);
    g_coset.numSims =numSims;

    % Randomizing targets
    targets = randomize_targets(g_coset);

    % Generating radar received analog signal (slow-fast matrix)
    x = generate_analog_input_signal(g_coset, targets);
    % x = stam(g, targets);

    % Processing
    targets_Coset = coset_nyquist(g_coset,x,2,targets);
    
%     if g.mf.debug_plot
%         for l=1:g.L
%             plot(targets.t(l), targets.f(l), 'bo');
%         end
%     end

    [isSuccess,realHist(i,:,:),resultHist(i,:,:),successVec] = analyze_result(g_coset,targets,targets_Coset,i);
    success = success + isSuccess;
    fprintf('Success rate: %.1f percent\n', 100*success/i);
   
end
if success == numSims
    fprintf('Great Success!\n');
end

% targets_Coset.f
% targets_Coset.t
% plot_results(g_coset, targets, targets_Coset);