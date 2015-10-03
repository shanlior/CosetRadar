function [successVec,resultHist,realHist,Results,targets,targets_Coset] = sim1()
%Ci= 30; % channel coefficient
%Q = 2; % How many ambiguities are resolved
% Ci = 79;
Q = 3;  
Ci=[3 7 8 11 91];
% Simulation config
rng('shuffle');
%r = rng(rngseed)
success = 0;
numSims = 10;
L = 2; % numTargets
realHist = zeros(numSims,L,2);
resultHist = zeros(numSims,L,2);
successVec = zeros(numSims,1);
for (i=1:numSims)

    Results(i).a=rng('shuffle');
    g = global_settings();
    g.P = 10;
    g_coset = global_settings_coset(g.P,g.P,L,-26, 1, 1, 0, 1,Ci,Q);

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

    targets_real = [round(targets.t /g_coset.CS.delta_t + 1) , round(targets.f *  g_coset.P * g_coset.tau + 1)];
    for l=1:L
        if targets_real(l,2) > g_coset.P
            targets_real(l,2) = targets_real(l,2) - 9;
        end
    end
    targets_real = sort(targets_real);
    targets_result = [targets_Coset.t , targets_Coset.f];
    targets_result = sort(targets_result);
    realHist(i,:,:) = targets_real(:,:);
    resultHist(i,:,:) = targets_result(:,:);
    %[targets_Coset,stats] = analyze_results(g_coset, targets, targets_Coset, 'NU_SubNyq');
    fprintf('iteration: %d\n', i);
    successVec(i) = max(max(abs(targets_real - targets_result)));
    
    if max(max(abs(targets_real - targets_result))) == 0
        success = success + 1;
    end
    fprintf('Success: %d\n', success);


end
success
% targets_Coset.f
% targets_Coset.t
%plot_results(g_coset, targets, targets_Coset);