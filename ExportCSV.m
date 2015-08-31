% headers = {'SNR', 'NU SubNyq FISTA Con. Coeff.', ...
%            'NU SubNyq 50/100 OMP Rand. Coeff.', ...
%            'NU SubNyq 20/100 OMP Rand. Coeff.', ...
%            'NU SubNyq 10/100 OMP Rand. Coeff.'};   
% 
% % Hit Rate
% % --------
% data = [SNR_dB_all(1:m)', hit_rate];
% csvwrite_with_headers('ResultsHitRateOMPRand.csv', data, headers)
% 
% % Time RMS Error
% % --------------
% data = [SNR_dB_all(1:m)', rms_t / g_NU_100.Nyquist.dt];
% csvwrite_with_headers('ResultsTimeRMSErrOMPRand.csv', data, headers)
% 
% % Frequency RMS Error
% % -------------------
% data = [SNR_dB_all(1:m)', rms_f / g_NU_100.Nyquist.df];
% csvwrite_with_headers('ResultsFreqRMSErrOMPRand.csv', data, headers)
% 
% % Time Error Bias
% % ---------------
% data = [SNR_dB_all(1:m)', bias_t / g_NU_100.Nyquist.dt];
% csvwrite_with_headers('ResultsTimeErrBiasOMPRand.csv', data, headers)
% 
% % Frequency Error Bias
% % ---------------------
% data = [SNR_dB_all(1:m)', bias_f / g_NU_100.Nyquist.df];
% csvwrite_with_headers('ResultsFreqErrBiasOMPRand.csv', data, headers)


% i = 1;
% dt_list = zeros(18,1);
% df_list = zeros(18,1);
% for m=100:5:190
%     g_tmp = global_settings_nu(100 ,m, 5, -26, 10, 1, 0, 1);
%     dt_list(i) = g_tmp.Nyquist.dt;
%     df_list(i) = g_tmp.Nyquist.df;
%     i=i+1;
% end
% 
% dt_list2 = [dt_list, dt_list];
% df_list2 = [df_list, df_list];
% 
% 
% headers = {'SNR', 'NU SubNyq OMP Con. Coeff.', ...
%            'NU SubNyq OMP Rand. Coeff.'};
%        
% % Hit Rate
% % --------
% data = [CPI', hit_rate];
% csvwrite_with_headers('ResultsHitRate_OMPCPISweep.csv', data, headers)
% 
% % Time RMS Error
% % --------------
% data = [CPI', rms_t ./ dt_list2];
% csvwrite_with_headers('ResultsTimeRMSErr_OMPCPISweep.csv', data, headers)
% 
% % Frequency RMS Error
% % -------------------
% data = [CPI', rms_f ./ df_list2];
% csvwrite_with_headers('ResultsFreqRMSErr_OMPCPISweep.csv', data, headers)
% 
% % Time Error Bias
% % ---------------
% data = [CPI', bias_t ./ dt_list2];
% csvwrite_with_headers('ResultsTimeErrBias_OMPCPISweep.csv', data, headers)
% 
% % Frequency Error Bias
% % ---------------------
% data = [CPI', bias_f ./ df_list2];
% csvwrite_with_headers('ResultsFreqErrBias_OMPCPISweep.csv', data, headers)

hit_rate_con = zeros(400,1);
for i=1:200
    hit_rate_con(4*i-3) = W_step(i).hit_rate.LPF_nu(1);
    hit_rate_con(4*i-2) = W_step(i).hit_rate.LPF_nu(2);
    hit_rate_con(4*i-1) = W_step(i).hit_rate.LPF_nu(3);
    hit_rate_con(4*i) = W_step(i).hit_rate.LPF_nu(4);
end

hit_rate_rand = zeros(400,1);
for i=1:200
    hit_rate_rand(4*i-3) = W_step(i).hit_rate.Random_nu(1);
    hit_rate_rand(4*i-2) = W_step(i).hit_rate.Random_nu(2);
    hit_rate_rand(4*i-1) = W_step(i).hit_rate.Random_nu(3);
    hit_rate_rand(4*i) = W_step(i).hit_rate.Random_nu(4);
end

csvwrite('Results_tmp_con5.csv', hit_rate_con);
csvwrite('Results_tmp_rand5.csv', hit_rate_rand);

% data = [];
% 
% for i=100:5:185
%     load(strcat('G:\CPI_', int2str(i)))
%     
%     data = [data; ones(15,1)*i, distances_f', P_separate_NU_SubNyq];
%     
%     clearvars -except data;
% end
% 
% csvwrite('CPISweep.csv', data);