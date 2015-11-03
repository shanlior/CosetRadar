function [success,realHist,resultHist,successVec] = analyze_result (g_coset,targets,targets_Coset,iterNum,plot_fail_sim)

success = 0;
successVec = zeros(1,2);


 targets_real = [round(targets.t /g_coset.CS.delta_t + 1) , round(targets.f *  g_coset.P * g_coset.tau + 1)];
    for l=1:g_coset.L
        if targets_real(l,2) > g_coset.P
            targets_real(l,2) = targets_real(l,2) - 9;
        end
    end
    targets_real_sort = sortrows(targets_real);
    targets_result = [targets_Coset.t , targets_Coset.f];
    targets_result_sort = sortrows(targets_result);
    realHist = targets_real(:,:);
    resultHist = targets_result(:,:);
    %[targets_Coset,stats] = analyze_results(g_coset, targets, targets_Coset, 'NU_SubNyq');
    differ=abs(targets_real_sort - targets_result_sort);
    successVec(1) = max(max(differ));
    for i = 1:g_coset.L 
        if differ(i,1) < 2 && differ(i,2) < 2 
            successVec(2) = successVec(2) + 1;
        end
    end
%     if successVec(1) == 9 
%         successVec(1) = 1;
%     end
    if successVec(1) < 2
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
            str_title{1}=sprintf('Failed : iter num %d , successVec(%d) = %d',iterNum,iterNum,successVec(1));
            title(str_title)
             ylabel('Time');
            xlabel('Frequency');
            legend('Real','Result')
            grid on
            hold off  
        end
    end
end
