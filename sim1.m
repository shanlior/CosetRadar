function [Results,targets,targets_Coset] = sim1()
%Ci= 30; % channel coefficient
%Q = 2; % How many ambiguities are resolved
% Ci = 79;
Q = 3;  
Ci=[79 83  89 97 101 103 113];
% Simulation config
rng(15);
%r = rng(rngseed)
success = 0;
for (i=1)

    Results(i).a=rng(17);
    g = global_settings();
    g_coset = global_settings_coset(g.P,100,2,-26, 1, 1, 0, 1,Ci,Q);

    % Randomizing targets
    targets = randomize_targets(g_coset);

    % Generating radar received analog signal (slow-fast matrix)
    x = generate_analog_input_signal(g_coset, targets);
    % x = stam(g, targets);

    % Processing
    targets_Coset = coset_nyquist(g_coset,x,2,targets);

    if g.mf.debug_plot
        for l=1:g.L
            plot(targets.t(l), targets.f(l), 'bo');
        end
    end

    % targets_Coset.f
    % targets_Coset.t
    %[targets_Coset,stats] = analyze_results(g_coset, targets, targets_Coset, 'NU_SubNyq');
%     disp iteration:
%     i
%     disp success:
%     Results(i).num_hits = stats.num_hits;
%     if (stats.num_hits == 5)
%         success = success + 1
%         Results(i).success = 1;
%     else
%         Results(i).success = 0;
%     end
end
% success
% targets_Coset.f
% targets_Coset.t
%plot_results(g_coset, targets, targets_Coset);