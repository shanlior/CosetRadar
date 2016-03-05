snr = -50:0;
numSims = 100;
L = 5;
for i=1:length(snr)
    [successVec,resultHist,realHist,targets,targets_Coset] = sim1([31 79],2,L,100,snr(i),false,numSims,1,1,1,100,100,true);
    hitrate(i) = sum(successVec)/(numSims * L);
    fprintf('Hitrate for snr = %d is %.2f\n',snr(i),hitrate(i));
end

save('hitrate.mat',hitrate);