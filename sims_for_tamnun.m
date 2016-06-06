close all;clear all;
L=5;
Q=2;
numSims = 100;
P=100;
snr = [-43 -40 -38 -35 -30 -25 -20];
snr = -50:0;
all_primes = primes(500);
all_primes = [0 all_primes];
save_opt=1;
config = 7;
%% config 1  :    Full Sample , num Ci = Q,  random Ci , all pulses - coset radar
if config == 1
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
                   sim1(Ci(j,:),Q,L,P,snr(i),false,numSims,1,1,1,100,100,0);
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
end



%% config 2  :    Full Sample , num Ci > Q,  rand Ci , partial  pulses and B , not same partial pulses and B- coset radar
if config == 2
    num_of_cfg = 4;
    partial_pulses = [50 25 20];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_2'];
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10};
      mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            for j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,1,1,1,partial_pulses(ii),partial_pulses(ii),0);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/Ci_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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


%% config 3  :    Full Sample , num Ci > Q,  rand Ci , partial  pulses and B , same partial pulses and B- coset radar
if config == 3
    num_of_cfg = 4;
    partial_pulses = [50 25 20];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_3'];
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    
    all_Ci = {Ci_4,Ci_8,Ci_10};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            for j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,1,1,1,partial_pulses(ii),partial_pulses(ii),1);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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

%% config 4  :    Full Sample , num Ci > Q,  rand Ci , partial  pulses , full B , not same partial pulses and B- coset radar
if config == 4
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_4'];
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,1,1,1,100,partial_pulses(ii),0);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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
    less_p = [75 50 25];
    success_rate_per_targets = zeros(length(less_p),length(snr));
    success_per_targets = zeros(length(less_p),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_6'];
    
    mkdir(cfg_path);
    for j = 1: length(less_p)
%         for i = 1:length(snr)
        parfor i = 1:length(snr)
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


%% config 8  :   not Full Sample = 25 , num Ci > Q,  rand Ci , partial  pulses and B , not same partial pulses and B- coset radar
if config == 8
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_8'];
    
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,0,1,4,partial_pulses(ii),partial_pulses(ii),0);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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


%% config 9  :   not Full Sample = 25  , num Ci > Q,  rand Ci , partial  pulses and B , same partial pulses and B- coset radar
if config == 9
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_9'];
    
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,0,1,4,partial_pulses(ii),partial_pulses(ii),1);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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

%% config 10  :    not Full Sample = 25 , num Ci > Q,  rand Ci , partial  pulses , full B , not same partial pulses and B- coset radar
if config == 10
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_10'];
    
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,0,1,4,100,partial_pulses(ii),0);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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


%% config 11  :   not Full Sample = 20 , num Ci > Q,  rand Ci , partial  pulses and B , not same partial pulses and B- coset radar
if config == 11
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_11'];
    
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,0,1,5,partial_pulses(ii),partial_pulses(ii),0);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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


%% config 12  :   not Full Sample = 20  , num Ci > Q,  rand Ci , partial  pulses and B , same partial pulses and B- coset radar
if config == 12
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_12'];
    
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,0,1,5,partial_pulses(ii),partial_pulses(ii),1);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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

%% config 13  :    not Full Sample = 20 , num Ci > Q,  rand Ci , partial  pulses , full B , not same partial pulses and B- coset radar
if config == 13
    num_of_cfg = 12;
    partial_pulses = [50 25 20 10];
    success_rate_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    success_per_targets = zeros(num_of_cfg,length(snr),size(partial_pulses,1));
    cur_dir=pwd;
    cfg_path = [cur_dir '/cfg_13'];
    
    load('Ci_4.mat');
    load('Ci_8.mat');
    load('Ci_10.mat');
    load('Ci_20.mat');
    all_Ci = {Ci_4,Ci_8,Ci_10,Ci_20};
    mkdir(cfg_path);
    for i = 1:length(snr)
        for ii = 1:length(partial_pulses)
            parfor j=1:num_of_cfg
                 tic
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                       sim1(all_Ci{ii}(j,:),Q,L,P,snr(i),false,numSims,0,1,5,100,partial_pulses(ii),0);
                success_per_targets(j,i,ii) = sum(successVec);
                success_rate_per_targets(j,i,ii) = 100*success_per_targets(j,i,ii) / size(successVec,1) / L; 
                if (save_opt) 
                    str_snr=int2str(snr(i));
                    str_j=int2str(j);
                    str_ii = int2str(partial_pulses(ii));
                    dest = [cfg_path '/iter_' str_j '/par_pulses' str_ii  '/snr_' str_snr '/'];
                    mkdir(dest);
                    f_dest = [dest 'successVec'];
                    parsave(f_dest,successVec);
                    f_dest = [dest 'success_per_targets'];
                    tmp_var = success_per_targets(j,i,ii);
                    parsave(f_dest,tmp_var);
                    f_dest = [dest 'success_rate_per_targets'];
                    tmp_var = success_rate_per_targets(j,i,ii);
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