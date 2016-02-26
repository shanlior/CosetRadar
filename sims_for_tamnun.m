close all;clear all;
L=5;
Q=2;
numSims = 100;
numSims = 1;
P=100;
snr = [-45 -40 -38 -35 -30 -25 -20];
all_primes = primes(500);
all_primes = [0 all_primes];
save_opt=1;
%% config 1  :    Full Sample , num Ci = Q,  random Ci  - coset radar
num_of_cfg = 24;
success_rate_per_targets = zeros(num_of_cfg,length(snr));
success_per_targets = zeros(num_of_cfg,length(snr));
Ci = zeros(num_of_cfg,Q);
for i=1:num_of_cfg
    Ci(i,:) = all_primes(randsample(length(all_primes),Q));
end
cur_dir=pwd;
cfg_1_path = [cur_dir '/cfg_1'];
mkdir(cfg_1_path);
save([cfg_1_path '/Ci'],'Ci');
for i = 1:length(snr)
    parfor j=1:num_of_cfg
         tic
         [successVec,resultHist,realHist,targets,targets_Coset] = ...
               sim1(Ci(j,:),Q,L,P,snr(i),false,numSims,1,1,1,100,100);
        success_per_targets(j,i) = sum(successVec);
        success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 
        if (save_opt) 
            str_snr=int2str(snr(i));
            str_j=int2str(j);
            dest = [cfg_1_path '/iter_' str_j '/snr_' str_snr '/'];
            mkdir(dest);
            f_dest = [dest 'successVec'];
            parsave(f_dest,successVec);
            f_dest = [dest 'success_per_targets'];
            tmp_var = success_per_targets(j,i);
            parsave(f_dest,tmp_var);
            f_dest = [dest 'success_rate_per_targets'];
            tmp_var = success_rate_per_targets(j,i);
            parsave(f_dest,tmp_var);
            f_dest = [dest 'resultHist'];
            parsave(f_dest,resultHist);
            f_dest = [dest 'realHist'];
            parsave(f_dest,realHist);
            f_dest = [dest 'targets'];
            parsave(f_dest,targets);
            f_dest = [dest 'targets_Coset'];
            parsave(f_dest,targets_Coset);
        end
        toc
    end
end
if save_opt 
    dest = [cfg_1_path '/'];
    f_dest = [dest 'success_per_targets'];
    tmp_var = success_per_targets;
    save(f_dest,'tmp_var');
    f_dest = [dest 'success_rate_per_targets'];
    tmp_var = success_rate_per_targets;
    save(f_dest,'tmp_var');
end

% figure
% plot(snr,success_rate_per_targets(1,:),'o-',snr,success_rate_per_targets(2,:),'x-')
% str_title{1}=sprintf('Coset: L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = detected targets / total targets');
% title(str_title,'FontSize',14)
% axis([snr(1)-1 snr(end)+1 -10 110]);
% xlabel('snr [dB]','FontSize',14);
% ylabel('Success Rate','FontSize',14);
% legend('Q=2','Q=3');
% set(gca,'FontSize',14);
