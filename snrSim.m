close all;clear all;
L=5;
numSims = 100;
P=100;
save_opt=1;
%% graph 1    success rate vs snr  - coset radar
% Q=[2,3];
% Ci=[0 3 5 7];
% snr=[-50:1:0];
% success_rate_per_sim_2 = zeros(length(Q),length(snr));
% success_rate_per_targets_2 = zeros(length(Q),length(snr));
% success_per_sim_2 = zeros(length(Q),length(snr));
% success_per_targets_2 = zeros(length(Q),length(snr));
% if save_opt
%     cur_dir=pwd;
%     main_result_path = [ cur_dir '\sims_results\coset'];
% end
% for j=1:length(Q)
%     for i = 1:length(snr)
%          fprintf('graph 1, snr = %d, Q=%d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!coset!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', snr(i),Q(j));
%          tic
%          [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci(1:Q(j)),Q(j),L,P,snr(i),false,numSims);
%          success_per_sim_2(j,i) = length(find(successVec==L));
%          success_rate_per_sim_2(j,i) = 100*success_per_sim_2(j,i) / size(successVec,1);
%         success_per_targets_2(j,i) = sum(successVec);
%         success_rate_per_targets_2(j,i) = 100*success_per_targets_2(j,i) / size(successVec,1) / L; 
%         if (save_opt) 
%             str_snr=int2str(snr(i));
%             str_j=int2str(j);
%             dest = [main_result_path '\Q_is_' str_j '\snr_is_' str_snr '\'];
%             mkdir(dest);
%             f_dest = [dest 'successVec'];
%             save(f_dest,'successVec');
%             f_dest = [dest 'success_per_sim'];
%             tmp_var = success_per_sim_2(j,i);
%             save(f_dest,'tmp_var');
%             f_dest = [dest 'success_rate_per_sim'];
%             tmp_var = success_rate_per_sim_2(j,i);
%             save(f_dest,'tmp_var');
%             f_dest = [dest 'success_per_targets'];
%             tmp_var = success_per_targets_2(j,i);
%             save(f_dest,'tmp_var');
%             f_dest = [dest 'success_rate_per_targets'];
%             tmp_var = success_rate_per_targets_2(j,i);
%             save(f_dest,'tmp_var');
%             f_dest = [dest 'resultHist'];
%             save(f_dest,'resultHist');
%             f_dest = [dest 'realHist'];
%             save(f_dest,'realHist');
%             f_dest = [dest 'targets'];
%             save(f_dest,'targets');
%             f_dest = [dest 'targets_Coset'];
%             save(f_dest,'targets_Coset');
%         end
%         toc
%     end
%     if save_opt
%         dest = [main_result_path '\Q_is_' str_j '\'];
%         f_dest = [dest 'success_per_sim'];
%         tmp_var = success_per_sim_2(j,:);
%         save(f_dest,'tmp_var');
%         f_dest = [dest 'success_rate_per_sim'];
%         tmp_var = success_rate_per_sim_2(j,:);
%         save(f_dest,'tmp_var');
%         f_dest = [dest 'success_per_targets'];
%         tmp_var = success_per_targets_2(j,:);
%         save(f_dest,'tmp_var');
%         f_dest = [dest 'success_rate_per_targets'];
%         tmp_var = success_rate_per_targets_2(j,:);
%         save(f_dest,'tmp_var');
%     end
% end
% if save_opt 
%     dest = [main_result_path '\'];
%     f_dest = [dest 'success_per_sim'];
%     tmp_var = success_per_sim_2;
%     save(f_dest,'tmp_var');
%     f_dest = [dest 'success_rate_per_sim'];
%     tmp_var = success_rate_per_sim_2;
%     save(f_dest,'tmp_var');
%     f_dest = [dest 'success_per_targets'];
%     tmp_var = success_per_targets_2;
%     save(f_dest,'tmp_var');
%     f_dest = [dest 'success_rate_per_targets'];
%     tmp_var = success_rate_per_targets_2;
%     save(f_dest,'tmp_var');
% end
% 
% figure
% plot(snr,success_rate_per_sim_2(1,:),'o-',snr,success_rate_per_sim_2(2,:),'x-')
% str_title{1}=sprintf('Coset: L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = sims passed / total sims');
% title(str_title,'FontSize',14)
% axis([snr(1)-1 snr(end)+1 -10 110]);
% xlabel('snr [dB]','FontSize',14);
% ylabel('Success Rate','FontSize',14);
% legend('Q=2','Q=3');
% set(gca,'FontSize',14);
% 
% figure
% plot(snr,success_rate_per_targets_2(1,:),'o-',snr,success_rate_per_targets_2(2,:),'x-')
% str_title{1}=sprintf('Coset: L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = detected targets / total targets');
% title(str_title,'FontSize',14)
% axis([snr(1)-1 snr(end)+1 -10 110]);
% xlabel('snr [dB]','FontSize',14);
% ylabel('Success Rate','FontSize',14);
% legend('Q=2','Q=3');
% set(gca,'FontSize',14);

%% success rate vs snr  - staggered radar
Q=[2,3];
Ci=[0 3 5 7];
snr=[-50:1:0];
success_rate_per_sim = zeros(length(Q),length(snr));
success_rate_per_targets = zeros(length(Q),length(snr));
success_per_sim = zeros(length(Q),length(snr));
success_per_targets = zeros(length(Q),length(snr));
if save_opt
    cur_dir=pwd;
    main_result_path = [ cur_dir '\sims_results\staggered'];
end
for j=1:length(Q)
    for i = 1:length(snr)
         fprintf('graph 1, snr = %d, Q=%d  !!!!!!!!!!!!!!!!!!!!!!!!!!staggered!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', snr(i),Q(j));
         tic
        [successVec,resultHist,realHist,targets,targets_staggered] = sim_staggered(0,Q(j),L,P,snr(i),false,numSims);
        success_per_sim(j,i) = length(find(successVec==L));
        success_rate_per_sim(j,i) = 100*success_per_sim(j,i) / size(successVec,1);
        success_per_targets(j,i) = sum(successVec);
        success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 
        if (save_opt) 
            str_snr=int2str(snr(i));
            str_j=int2str(j);
            dest = [main_result_path '\Q_is_' str_j '\snr_is_' str_snr '\'];
            mkdir(dest);
            f_dest = [dest 'successVec'];
            save(f_dest,'successVec');
            f_dest = [dest 'success_per_sim'];
            tmp_var = success_per_sim(j,i);
            save(f_dest,'tmp_var');
            f_dest = [dest 'success_rate_per_sim'];
            tmp_var = success_rate_per_sim(j,i);
            save(f_dest,'tmp_var');
            f_dest = [dest 'success_per_targets'];
            tmp_var = success_per_targets(j,i);
            save(f_dest,'tmp_var');
            f_dest = [dest 'success_rate_per_targets'];
            tmp_var = success_rate_per_targets(j,i);
            save(f_dest,'tmp_var');
            f_dest = [dest 'resultHist'];
            save(f_dest,'resultHist');
            f_dest = [dest 'realHist'];
            save(f_dest,'realHist');
            f_dest = [dest 'targets'];
            save(f_dest,'targets');
            f_dest = [dest 'targets_staggered'];
            save(f_dest,'targets_staggered');
        end
        toc
    end
    if save_opt
        dest = [main_result_path '\Q_is_' str_j '\'];
        f_dest = [dest 'success_per_sim'];
        tmp_var = success_per_sim(j,:);
        save(f_dest,'tmp_var');
        f_dest = [dest 'success_rate_per_sim'];
        tmp_var = success_rate_per_sim(j,:);
        save(f_dest,'tmp_var');
        f_dest = [dest 'success_per_targets'];
        tmp_var = success_per_targets(j,:);
        save(f_dest,'tmp_var');
        f_dest = [dest 'success_rate_per_targets'];
        tmp_var = success_rate_per_targets(j,:);
        save(f_dest,'tmp_var');
    end
end
if save_opt 
    dest = [main_result_path '\'];
    f_dest = [dest 'success_per_sim'];
    tmp_var = success_per_sim;
    save(f_dest,'tmp_var');
    f_dest = [dest 'success_rate_per_sim'];
    tmp_var = success_rate_per_sim;
    save(f_dest,'tmp_var');
    f_dest = [dest 'success_per_targets'];
    tmp_var = success_per_targets;
    save(f_dest,'tmp_var');
    f_dest = [dest 'success_rate_per_targets'];
    tmp_var = success_rate_per_targets;
    save(f_dest,'tmp_var');
end

figure
plot(snr,success_rate_per_sim(1,:),'o-',snr,success_rate_per_sim(2,:),'x-')
str_title{1}=sprintf('Staggered: L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = sims passed / total sims');
title(str_title,'FontSize',14)
axis([snr(1)-1 snr(end)+1 -10 110]);
xlabel('snr [dB]','FontSize',14);
ylabel('Success Rate','FontSize',14);
legend('Q=2','Q=3');
set(gca,'FontSize',14);

figure
plot(snr,success_rate_per_targets(1,:),'o-',snr,success_rate_per_targets(2,:),'x-')
str_title{1}=sprintf('Staggered: L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = detected targets / total targets');
title(str_title,'FontSize',14)
axis([snr(1)-1 snr(end)+1 -10 110]);
xlabel('snr [dB]','FontSize',14);
ylabel('Success Rate','FontSize',14);
legend('Q=2','Q=3');
set(gca,'FontSize',14);
