close all;clear all;
L=5;
numSims = 2;
P=100;
%% graph 1    success rate vs snr 
Q=[1,2,4];
Ci=[0 3 5 7];
snr=[-50:2.5:0];
success_rate_per_sim_2 = zeros(length(Q),length(snr));
success_rate_per_targets_2 = zeros(length(Q),length(snr));
success_per_sim_2 = zeros(length(Q),length(snr));
success_per_targets_2 = zeros(length(Q),length(snr));
for j=1:length(Q)
    for i = 1:length(snr)
         fprintf('graph 1, snr = %d, Q=%d  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', snr(i),Q(j));
         [successVec,resultHist,realHist,targets,targets_Coset] = sim1(Ci(1:Q(j)),Q(j),L,P,snr(i),false,numSims);
         success_per_sim_2(j,i) = length(find(successVec(:,1) < 2));
         success_rate_per_sim_2(j,i) = 100*success_per_sim_2(j,i) / size(successVec,1);
        success_per_targets_2(j,i) = sum(successVec(:,2));
        success_rate_per_targets_2(j,i) = 100*success_per_targets_2(j,i) / size(successVec,1) / L; 
    end
end

figure
plot(snr,success_rate_per_sim_2(1,:),'o-',snr,success_rate_per_sim_2(2,:),'x-')
str_title{1}=sprintf('L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = sims passed / total sims');
title(str_title,'FontSize',14)
axis([snr(1)-1 snr(end)+1 -10 110]);
xlabel('snr [dB]','FontSize',14);
ylabel('Success Rate','FontSize',14);
legend('Q=2','Q=4');
set(gca,'FontSize',14);

figure
plot(snr,success_rate_per_targets_2(1,:),'o-',snr,success_rate_per_targets_2(2,:),'x-')
str_title{1}=sprintf('L = 5 , P = 100 ,length(Ci)=Q , Success rate vs Noise\n rate = detected targets / total targets');
title(str_title,'FontSize',14)
axis([snr(1)-1 snr(end)+1 -10 110]);
xlabel('snr [dB]','FontSize',14);
ylabel('Success Rate','FontSize',14);
legend('Q=2','Q=4');
set(gca,'FontSize',14);


