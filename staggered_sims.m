%%sim staggered

config = 10;

save_opt
T = 2;
L = 5;
P = 100;
Q = 2;
snr = -50:0;
sample_SubNyquist_factor = 1;
pulse_SubNyquist_factor = 1;


%% config 1  :   
if config == 1 || config == 10
    cfg5_Q = [2 4];
    sample_SubNyquist_factor = 1;
    pulse_SubNyquist_factor = 1;
    focusing = 0;
    success_rate_per_targets = zeros(length(cfg5_Q),length(snr));
    success_per_targets = zeros(length(cfg5_Q),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/staggered_less_pulses'];
    
    mkdir(cfg_path);
    for j = 1:length(cfg5_Q)
        parfor i = 1:length(snr)
            str_line = ['-------cfg5-------Q = ',num2str(cfg5_Q(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_staggered] = sim_staggered(0,Q,L,P,P,snr(i),false,100,sample_SubNyquist_factor, cfg5_Q(j),focusing);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 

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
    
end


%% config 2  :   less samples
if config == 2 || config == 10
    cfg2_Q = [2 4 8];
    sample_SubNyquist_factor = 1;
    pulse_SubNyquist_factor = 1;
    focusing = 0;
    success_rate_per_targets = zeros(length(cfg2_Q),length(snr));
    success_per_targets = zeros(length(cfg2_Q),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/staggered_less_samples'];
    
    mkdir(cfg_path);
    for j = 1:length(cfg2_Q)
        parfor i = 1:length(snr)
            str_line = ['-------cfg2-------Q = ',num2str(cfg2_Q(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_staggered] = sim_staggered(0,Q,L,P,P,snr(i),false,100,cfg2_Q(j),pulse_SubNyquist_factor,focusing);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 

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
    
end

%% config 3  :   
if config == 3 || config == 9
    cfg5_Q = [2 4];
    sample_SubNyquist_factor = 1;
    pulse_SubNyquist_factor = 1;
    focusing = 1;
    success_rate_per_targets = zeros(length(cfg5_Q),length(snr));
    success_per_targets = zeros(length(cfg5_Q),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/staggered_less_pulses_focusing'];
    
    mkdir(cfg_path);
    for j = 1:length(cfg5_Q)
        parfor i = 1:length(snr)
            str_line = ['-------cfg5-------Q = ',num2str(cfg5_Q(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_staggered] = sim_staggered(0,Q,L,P,P,snr(i),false,100,sample_SubNyquist_factor, cfg5_Q(j),focusing);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 

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
    
end


%% config 4  :   less samples
if config == 4 || config == 10
    cfg2_Q = [2 4 8];
    sample_SubNyquist_factor = 1;
    pulse_SubNyquist_factor = 1;
    focusing = 1;
    success_rate_per_targets = zeros(length(cfg2_Q),length(snr));
    success_per_targets = zeros(length(cfg2_Q),length(snr));
    cur_dir=pwd;
    cfg_path = [cur_dir '/staggered_less_samples_focusing'];
    
    mkdir(cfg_path);
    for j = 1:length(cfg2_Q)
        parfor i = 1:length(snr)
            str_line = ['-------cfg2-------Q = ',num2str(cfg2_Q(j)),'-------------snr = ',num2str(snr(i)),'------------'];
            disp(str_line)
             tic
             [successVec,resultHist,realHist,targets,targets_staggered] = sim_staggered(0,Q,L,P,P,snr(i),false,100,cfg2_Q(j),pulse_SubNyquist_factor,focusing);
            success_per_targets(j,i) = sum(successVec);
            success_rate_per_targets(j,i) = 100*success_per_targets(j,i) / size(successVec,1) / L; 

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
    
end