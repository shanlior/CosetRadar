close all;clear all;
L=5;
Q=2;
numSims = 1;
P=100;
snr = [-43 -40 -38 -35 -30 -25 -20];
all_primes = primes(500);
all_primes = [0 all_primes];
save_opt=1;
config = 2;
matlabpool(12)
tic
parfor i=1:12
             [successVec,resultHist,realHist,targets,targets_Coset] = ...
                   sim1([2 3],2,5,100,-30,false,numSims,1,1,1,100,100,0);
               
end
disp 'Parallel loop'
toc
matlabpool close
tic
for i=1:12
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                   sim1([2 3],2,5,100,-30,false,numSims,1,1,1,100,100,0);
end
disp 'Serial loop'
toc
matlabpool(12)

tic
for i=1:12
                 [successVec,resultHist,realHist,targets,targets_Coset] = ...
                   sim1([2 3],2,5,100,-30,false,numSims,1,1,1,100,100,0);
end
disp 'Serial Pool loop'
toc
matlabpool close

