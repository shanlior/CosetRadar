close all;clear all;
%% graph 1
P=10;
Ci=[0 3 5 7 11  17 19 23 ];
max_num_targets = 4;
success_rate_per_sim_1 = zeros(1,max_num_targets);
success_rate_per_targets_1 = zeros(1,max_num_targets);
success_per_sim_1 = zeros(1,max_num_targets);
success_per_targets_1 = zeros(1,max_num_targets);
for Q = [2,4]
    for L = 1:max_num_targets
        fprintf('graph 1, Q = %d ; L = %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', Q,L);
        [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,false,100);
        success_per_sim_1(L) = length(find(successVec(:,1) < 2));
        success_rate_per_sim_1(L) = 100*success_per_sim_1(L) / length (successVec);
        success_per_targets_1(L) = sum(successVec(:,2));
        success_rate_per_targets_1(L) = 100*success_per_targets_1(L) / length (successVec) / L;      
    end
    
    figure
    subplot(2,1,1); plot(1:max_num_targets,success_rate_per_sim_1,'o-')
    str_title{1}=sprintf('Ci=[0 3 5 7 11  17 19 23 ]  ;  Q = %d ; P = 10 \n rate = sims passed / total sims',Q,'FontSize',14);
    title(str_title,'FontSize',14,'FontWeight','bold')
    axis([0.9 max_num_targets+0.1 -10 110]);
    xlabel('Number of Targets - L','FontSize',14);
    ylabel('Success Rate','FontSize',14);
    ax = gca;
    set(ax , 'XTick',1:max_num_targets,'XTickLabel',);
    subplot(2,1,2); plot(1:max_num_targets,success_rate_per_targets_1,'o-')
    str_title{1}=sprintf('rate = detected targets / total targets');
    title(str_title,'FontSize',14,'FontWeight','bold')
    axis([0.9 max_num_targets+0.1 -10 110]);
    ax = gca;
    set(ax , 'XTick',1:max_num_targets)
    xlabel('Number of Targets - L','FontSize',14);
    ylabel('Success Rate','FontSize',14);
end

%% graph 2
Q=2;
L=2;
P=10;
Ci=[0 3 5 7 11  17 19 23 ];
success_rate_per_sim_2 = zeros(1,length(Ci));
success_rate_per_targets_2 = zeros(1,length(Ci));
success_per_sim_2 = zeros(1,length(Ci));
success_per_targets_2 = zeros(1,length(Ci));

for i = 1:length(Ci)
     fprintf('graph 2, Ci length = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', i);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci(1:i),Q,L,P,false,100);
     success_per_sim_2(i) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_2(i) = 100*success_per_sim_2(i) / length (successVec);
    success_per_targets_2(i) = sum(successVec(:,2));
    success_rate_per_targets_2(i) = 100*success_per_targets_2(i) / length (successVec) / L; 
end

figure
subplot(2,1,1); plot(1:length(Ci),success_rate_per_sim_2,'o-')
str_title{1}=sprintf('L = 2 ; Q = 2 , P = 10 \n rate = sims passed / total sims');
title(str_title)
axis([0.9 length(Ci)+0.1 -10 110]);
xlabel('Number of Channels');
ylabel('Success Rate');
ax = gca;
set(ax , 'XTick',1:length(Ci));
subplot(2,1,2); plot(1:length(Ci),success_rate_per_targets_2,'o-')
str_title{1}=sprintf('rate = detected targets / total targets');
title(str_title)
axis([0.9 length(Ci)+0.1 -10 110]);
ax = gca;
set(ax , 'XTick',1:length(Ci))
xlabel('Number of Channels');
ylabel('Success Rate');


%% graph 3
P=10;
L=4;
Ci=[0 3 5 7 11  13 17 19 ];
success_rate_per_sim_3 = zeros(1,4);
success_rate_per_targets_3 = zeros(1,4);
success_per_sim_3 = zeros(1,4);
success_per_targets_3 = zeros(1,4);
for Q = 1:4
     fprintf('graph 3,  Q = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', Q);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,false,100);
     success_per_sim_3(Q) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_3(Q) = 100*success_per_sim_3(i) / length (successVec);
    success_per_targets_3(Q) = sum(successVec(:,2));
    success_rate_per_targets_3(Q) = 100*success_per_targets_3(i) / length (successVec) / L; 
end
figure
subplot(2,1,1); plot(1:4,success_rate_per_sim_3,'o-')
str_title{1}=sprintf('L = 4 ; Ci = [0 3 5 7 11 13 17 19] ; P = 10 \n rate = sims passed / total sims');
title(str_title)
axis([0.9 4+0.1 -10 110]);
xlabel('Q');
ylabel('Success Rate');
ax = gca;
set(ax , 'XTick',1:4);
subplot(2,1,2); plot(1:4,success_rate_per_targets_3,'o-')
str_title{1}=sprintf('rate = detected targets / total targets');
title(str_title)
axis([0.9 4+0.1 -10 110]);
ax = gca;
set(ax , 'XTick',1:4)
xlabel('Q');
ylabel('Success Rate');

%% graph 4
Q=4;
L=2;
P=10;
Ci=[0 3 5 7 11  17 19 23 ];
success_rate_per_sim_4 = zeros(1,length(Ci));
success_rate_per_targets_4 = zeros(1,length(Ci));
success_per_sim_4 = zeros(1,length(Ci));
success_per_targets_4 = zeros(1,length(Ci));

for i = 1:length(Ci)
     fprintf('graph 2, Ci length = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', i);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci(1:i),Q,L,P,false,100);
     success_per_sim_4(i) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_4(i) = 100*success_per_sim_4(i) / length (successVec);
    success_per_targets_4(i) = sum(successVec(:,2));
    success_rate_per_targets_4(i) = 100*success_per_targets_4(i) / length (successVec) / L; 
end

figure
subplot(2,1,1); plot(1:length(Ci),success_rate_per_sim_4,'o-')
str_title{1}=sprintf('L = 2 ; Q = 4 , P = 10 \n rate = sims passed / total sims');
title(str_title)
axis([0.9 length(Ci)+0.1 -10 110]);
xlabel('Number of Channels');
ylabel('Success Rate');
ax = gca;
set(ax , 'XTick',1:length(Ci));
subplot(2,1,2); plot(1:length(Ci),success_rate_per_targets_4,'o-')
str_title{1}=sprintf('rate = detected targets / total targets');
title(str_title)
axis([0.9 length(Ci)+0.1 -10 110]);
ax = gca;
set(ax , 'XTick',1:length(Ci))
xlabel('Number of Channels');
ylabel('Success Rate');

%% graph 5 part A
Q=3;
L=4;
P=10;
snr_db = [-20 -15 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2  -1 0 5 10 20];
Ci=[0 3 5 7 11];
success_rate_per_sim_5 = zeros(1,length(Ci));
success_rate_per_targets_5 = zeros(1,length(Ci));
success_per_sim_5 = zeros(1,length(Ci));
success_per_targets_5 = zeros(1,length(Ci));

for i = 1:length(snr_db)
     fprintf('graph 2, snr_db_number = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', i);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,snr_db(i),false,20);
     success_per_sim_5(i) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_5(i) = 100*success_per_sim_5(i) / length (successVec);
    success_per_targets_5(i) = sum(successVec(:,2));
    success_rate_per_targets_5(i) = success_per_targets_5(i) / length (successVec) / L; 
end

figure
plot(snr_db,success_rate_per_targets_5,'o-')
str_title{1}=sprintf('L = 4 ; Q = 2 , P = 10 \n rate = sims passed / total sims');
title(str_title)
axis([-20 20 0 1]);
xlabel('SNR');
ylabel('Hit Rate');
% ax = gca;
% set(ax , 'XTick',1:length(Ci));

Q=1;
L=4;
P=10;
snr_db = [-20 -15 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2  -1 0 5 10 20];
Ci=[0 3 5 7 11  17 19 23 ];
success_rate_per_sim_6 = zeros(1,length(Ci));
success_rate_per_targets_6 = zeros(1,length(Ci));
success_per_sim_6 = zeros(1,length(Ci));
success_per_targets_6 = zeros(1,length(Ci));

for i = 1:length(snr_db)
     fprintf('graph 2, snr_db_number = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', i);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,snr_db(i),false,20);
     success_per_sim_6(i) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_6(i) = 100*success_per_sim_6(i) / length (successVec);
    success_per_targets_6(i) = sum(successVec(:,2));
    success_rate_per_targets_6(i) = success_per_targets_6(i) / length (successVec) / L; 
end

figure
plot(snr_db,success_rate_per_targets_6,'o-')
str_title{1}=sprintf('L = 4 ; Q = 4 , P = 10 \n rate = sims passed / total sims');
title(str_title)
axis([-20 20 0 1]);
xlabel('SNR');
ylabel('Hit Rate');
% ax = gca;
% set(ax , 'XTick',1:length(Ci));





%% graph 5 part B
Q=3;
L=4;
P=10;
snr_db = [-20 -15 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2  -1 0 5 10 20];
Ci=[0 3 5 7 11];
success_rate_per_sim_7 = zeros(1,length(Ci));
success_rate_per_targets_7 = zeros(1,length(Ci));
success_per_sim_7 = zeros(1,length(Ci));
success_per_targets_7 = zeros(1,length(Ci));

for i = 1:length(snr_db)
     fprintf('graph 2, snr_db_number = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', i);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,snr_db(i),false,20);
     success_per_sim_7(i) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_7(i) = 100*success_per_sim_7(i) / length (successVec);
    success_per_targets_7(i) = sum(successVec(:,2));
    success_rate_per_targets_7(i) = success_per_targets_7(i) / length (successVec) / L; 
end

figure
plot(snr_db,success_rate_per_targets_7,'o-')
str_title{1}=sprintf('L = 4 ; Q = 3 , P = 10 \n Hit Rate vs. SNR');
title(str_title)
axis([-20 20 0 1]);
xlabel('SNR');
ylabel('Hit Rate');
% ax = gca;
% set(ax , 'XTick',1:length(Ci));

Q=1;
L=4;
P=10;
snr_db = [-20 -15 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2  -1 0 5 10 20];
Ci=[0 3 5 7 11  17 19 23 ];
success_rate_per_sim_8 = zeros(1,length(Ci));
success_rate_per_targets_8 = zeros(1,length(Ci));
success_per_sim_8 = zeros(1,length(Ci));
success_per_targets_8 = zeros(1,length(Ci));

for i = 1:length(snr_db)
     fprintf('graph 2, snr_db_number = %d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', i);
     [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci,Q,L,P,snr_db(i),false,20);
     success_per_sim_8(i) = length(find(successVec(:,1) < 2));
     success_rate_per_sim_8(i) = 100*success_per_sim_8(i) / length (successVec);
    success_per_targets_8(i) = sum(successVec(:,2));
    success_rate_per_targets_8(i) = success_per_targets_8(i) / length (successVec) / L; 
end

figure
plot(snr_db,success_rate_per_targets_8,'o-')
str_title{1}=sprintf('L = 4 ; Q = 1 , P = 10 \n Hit Rate vs. SNR');
title(str_title)
axis([-20 20 0 1]);
xlabel('SNR');
ylabel('Hit Rate');
% ax = gca;
% set(ax , 'XTick',1:length(Ci));

figure

plot(snr_db,success_rate_per_targets_5,'+-',snr_db,success_rate_per_targets_6,'o-',snr_db,success_rate_per_targets_7,'x-',snr_db,success_rate_per_targets_8,'*-')
title('Hit Rate vs. SNR');
axis([-20 20 0 1]);
xlabel('SNR');
ylabel('Hit Rate');
legend('Q=2','Q=4','Q=3','Q=1');



%% graph 6


[successVec,resultHist,realHist,targets,targets_Coset] = sim1();
