function [success,realHist,resultHist,successVec] = analyze_result (g_coset,targets,targets_Coset,iterNum)

success = 0;


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
    successVec(iterNum) = max(max(abs(targets_real_sort - targets_result_sort)));
    if successVec(iterNum) == 9 
        successVec(iterNum) = 1;
    end
    if successVec(iterNum) < 2
        success = success + 1;
        fprintf('Iteration %d is successful\n', iterNum);
    else
        fprintf('Iteration %d failed\n', iterNum);
    %%plot results Vs real   when failing
        figure
        hold on
        scatter(resultHist(:,2),resultHist(:,1))
        scatter(realHist(:,2),realHist(:,1),'.')
        axis([0 g_coset.P 0 g_coset.CS.N_t*g_coset.Q])
        str_title{1}=sprintf('Failed : iter num %d , successVec(%d) = %d',iterNum,iterNum,successVec(iterNum));
        title(str_title)
         ylabel('Time');
        xlabel('Frequency');
        legend('Result','Real')
        grid on
        hold off  
    end
end
