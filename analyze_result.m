function [success,realHist,resultHist,successVec] = analyze_result (g_coset,targets,targets_Coset,iterNum,plot_fail_sim,timeFlag)


if nargin~=6
    timeFlag = false;
end

if timeFlag
    targets_Coset.t = round(targets_Coset.t / g_coset.CS.delta_t + 1);
    targets_Coset.f = round(targets_Coset.f *  g_coset.P * g_coset.tau + 1);
end
success = 0;


 targets_real = [round(targets.t /g_coset.CS.delta_t + 1) , round(targets.f *  g_coset.P * g_coset.tau + 1)];
 
%     for l=1:g_coset.L
%         if targets_real(l,2) > g_coset.P
%             targets_real(l,2) = targets_real(l,2) - 9;
%         end
%     end


    P = perms(1:g_coset.L);
    hits = zeros(size(P));
    goodness = zeros(size(P,1),1);
%     R = diag([g_coset.hit_rate_threshold.t^-2 g_coset.hit_rate_threshold.f^-2]);
    R = diag([round(g_coset.hit_rate_threshold.t/g_coset.CS.delta_t)^-2 * round(g_coset.hit_rate_threshold.f * g_coset.P * g_coset.tau)^-2]);
%     R = eye(2);
    for p=1:size(P,1)
        t_error = targets_real(:,1) - targets_Coset.t(P(p,:));
        f_error = targets_real(:,2) - targets_Coset.f(P(p,:));

        d = [t_error f_error].';
        D = diag(d'*R*d);
        hits(p,:) = D'<1;
        goodness(p) = sum(1./(D(D<1)));
    end
    best_perm_candidates = find(sum(hits,2) == max(sum(hits,2)));
    [~, best_perm_index] = max(goodness(best_perm_candidates));
    best_perm_index = best_perm_candidates(best_perm_index);
    targets_Coset.best_perm = P(best_perm_index, :);
    targets_Coset.best_perm_hits = hits(best_perm_index, :) > 0;
    num_hits = sum(targets_Coset.best_perm_hits);
    
    successVec = num_hits;

% 
% 
%     targets_real_sort = sortrows(targets_real);
    targets_result = [targets_Coset.t , targets_Coset.f];
%     targets_result_sort = sortrows(targets_result);
    realHist = targets_real(:,:);
    resultHist = targets_result(:,:);
    %[targets_Coset,stats] = analyze_results(g_coset, targets, targets_Coset, 'NU_SubNyq');
%     differ=abs(targets_real_sort - targets_result_sort);
%     successVec(1) = max(max(differ));
%     for i = 1:g_coset.L 
%         if differ(i,1) < 2 && differ(i,2) < 2 
%             successVec(2) = successVec(2) + 1;
%         end
%     end
%     if successVec(1) == 9 
%         successVec(1) = 1;
%     end
    if successVec == g_coset.L
        success = success + 1;
        fprintf('Iteration %d is successful\n', iterNum);
    else
        fprintf('Iteration %d failed\n', iterNum);
    %%plot results Vs real   when failing
        if plot_fail_sim 
            figure
            hold on
            scatter(realHist(:,2),realHist(:,1))
            scatter(resultHist(:,2),resultHist(:,1),'.')
            axis([0 g_coset.P 0 g_coset.CS.N_t*g_coset.Q])
            str_title{1}=sprintf('Failed : iter num %d , successVec(%d) = %d',iterNum,iterNum,successVec);
            title(str_title)
             ylabel('Time');
            xlabel('Frequency');
            legend('Real','Result')
            grid on
            hold off  
        end
    end
end
