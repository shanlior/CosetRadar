close all;clear all;
L=5;
Q=2;
numSims = 100;
P=100;
snr = -50:0;
all_primes = primes(500);
all_primes = [0 all_primes];
save_opt=1;
config = 6;
%% config 5  :   Q = [3 4 8] ; Ci = [37 79] ; 
if config == 5
%     cfg5_Q = [2 3 4 8];
    cfg5_Q = [2 3 4 8];
    success_rate_per_targets = zeros(length(cfg5_Q),length(snr));
    success_per_targets = zeros(length(cfg5_Q),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_5'];
    
    mkdir(cfg_path);
    for j = 1: length(cfg5_Q)
        for i = 1:length(snr)
            str_line = ['-------cfg5-------Q = ',num2str(cfg5_Q(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_Coset] = ...
                   sim1([37 79],cfg5_Q(j),L,P,snr(i),false,numSims,1,1,1,P,P,1);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 
%             successVec = 0;
%             resultHist = 0;
%             realHist = 0;
%             targets = 0;
%             targets_Coset = 0;
            toc
        end
    end
    if save_opt 
        dest = [cfg_path '/'];
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
end



%% config 6  :   Q = 2 ; Ci = [37 79] ; less pulses = [25 50 75]
if config == 6
%     less_p = [75 50 25];
    less_p = [25];
    success_rate_per_targets = zeros(length(less_p),length(snr));
    success_per_targets = zeros(length(less_p),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_6'];
    
    mkdir(cfg_path);
    for j = 1: length(less_p)
        for i = 1:length(snr)
%         parfor i = 1:length(snr)
            str_line = ['-------cfg6-------lees_p = ',num2str(less_p(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_Coset] = ...
                   sim1([37 79],Q,L,P,snr(i),false,numSims,1,1,1,P,less_p(j),0);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 
%             successVec = 0;
%             resultHist = 0;
%             realHist = 0;
%             targets = 0;
%             targets_Coset = 0;
            toc
        end
    end
    if save_opt 
        dest = [cfg_path '/'];
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
end



%% config 7  :    not Full Sample
if config == 7
    sample_subNyquist_factor = [2 4 8];
    success_rate_per_targets = zeros(length(sample_subNyquist_factor),length(snr));
    success_per_targets = zeros(length(sample_subNyquist_factor),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/CosetRadarSybNyquistSamples'];
    
    mkdir(cfg_path);
    for j = 1: length(sample_subNyquist_factor)
%         for i = 1:length(snr)
        parfor i = 1:length(snr)
            str_line = ['-------cfg6-------Samples = ',num2str(sample_subNyquist_factor(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_Coset] = ...
                   sim1([37 79],Q,L,P,snr(i),false,numSims,0,1,sample_subNyquist_factor(j),P,P,0);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 
%             successVec = 0;
%             resultHist = 0;
%             realHist = 0;
%             targets = 0;
%             targets_Coset = 0;
            toc
        end
    end
    if save_opt 
        dest = [cfg_path '/'];
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
end
