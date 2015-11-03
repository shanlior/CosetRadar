
function [X, R, X_SR] = SolveOmpSingleLayerPhaseShiftFix2(Y, A1, B,A2 ,B2, numIters,sSimParams,sTargetsParams)
% solves Y = A*X*B.'
% Input:
%   Y, A, B     - matrices of equation
%   numIters    - number of algorithm iterations
% Output:
%   sResultsOmp - results by OMP algorithm


% 1. Initialization
% ------------------
R = Y; % < Residual
Supp = [];  % Support
SuppMat=[];
t = 1;
%A1_tag = A1';
%B_tag = B';
%A2_tag = A2';
[~,N] = size(A1);
%[N2,~] = size(A2);
[M,~] = size(B);
sumGetRid=0;


while t<=numIters %| norm(vec(R))> 80
%     str = ['Iteration #', num2str(t)];
%     disp(str);
    
    % 2. Find the two indicies of the support
    % ---------------------------------------
    %ThetaBin = -1 + sSimParams.sResolution.sinThetaBin*250;  %* (randi(MOrig * NOrig, L, 1)-1)
    
    % f=A_tag*R;
    %tic
%     for iChannels=1:M
%        %vecTheta=diag(B_tag(:,iChannels));
%        tt=B_tag(:,iChannels);
%        xx=length(tt);
%        tmp1 = abs(A_tag*(R*(sparse(1:xx,1:xx,tt)*S_tag))).';
%        tmp(:,iChannels)=vec(tmp1);
%     end
    %tic
    if sSimParams.fixAppertreDelays
        tmp=(A2'.*(A1'*R))*(B.*B2)';
    else
        tmp=(A2'.*(A1'*R))*B';
    end
    %toc
    
    %debug mode
    %stem3(abs(tmp))
   %view(180,0)
    
    
% stem3(abs(tmp))
% view(180,0)
% tmp1=A1'*R;
% stem3(abs(tmp1))
% view(180,0)
% view(90,0)
% stem3(abs(tmp1(1:1500,:)))
% view(90,0)
% stem3(abs(tmp1(1400:1600,:)))
% view(90,0)
% tmp2=A1'*R*B';
% stem3(abs(tmp2))
% view(90,0)
  
    %tmp1 = abs(R*B_tag(:,250));
       [~,i_vec] = max(tmp);
       [max_val, j] = max(max(tmp));
       i = i_vec(j);
    
    % 3. Augment index set
    % --------------------
   Supp(t,2) = j;  %sin
   Supp(t,1) = i;  %range 
   estimatedTau= (i-1)*sSimParams.sResolution.timeBinSuperResolution;
   estimatedTheta=-1+2*(j-1)/M;
   
%    j
%    mod(i,16)
%    i/16
%    www=1;
    
    
%     Supp_uniq = unique(Supp, 'rows');
%     while(size(Supp_uniq,1)<size(Supp,1))
%         %%disp('Support error');
%         tmp1(i,j)=0;
%         
%         [~,i_vec] = max(tmp1);
%         [max_val, j] = max(max(tmp1));
%         i = i_vec(j);  
%         
%         Supp(t,2) = j;
%         Supp(t,1) = i;
%         
%         Supp_uniq = unique(Supp, 'rows');
%     end
    
    % 4. Find the new signal estimate
    % -------------------------------
    %calculate matrix for the new entry

%     kVec                           = sSimParams.kappaSingleTrans; %+ ceil(sSimParams.NTag/2);
%     kVec                           = kVec(:);
%     freqShiftRatio                 = sSimParams.freqShiftRatio;
%     singleTransmitterBw            = sSimParams.singleTransmitterBw;
%     estimatedTau= (i-1)*sSimParams.sResolution.timeBinSuperResolution-sSimParams.sResolution.timeBin/2
%     %estimatedTau= 1.0e-06 * 0.4875
%     estimatedTauVec          = exp(-1j * 2 * pi / sSimParams.T * kVec * estimatedTau);
%     estimatedSRVec= exp(-1j * 2 * pi *freqShiftRatio * sSimParams.fcFreqDivisionLocation*estimatedTau);%
%     estimatedSRVec=repmat(estimatedSRVec,1,sSimParams.N);
%     aa=sSimParams.sResolution.timeBinSuperResolution;
%     i;
% 
%     
%     
%     estimatedTheta=-1+2*(j-1)/M
%     estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:));
%     
   
    
%     rr=S(mod(i-1,N2)+1,:);
%     yAxis=tmp((i-1):(i+1),j);
%     xAxis= [-1 0 1];
%     val = estimatedTau; %value to find
%     tmp1 = abs(sTargetsParams.tau-val);
%     [val2 idx] = min(tmp1); %index of closest value
%     closest = sTargetsParams.tau(idx); %closest value
    
%     C = [xAxis.^2;xAxis;ones(1,3)]\yAxis;
%     
%     midValEstParbolicFit= -C(2)/(2*C(1));
%     
%     C(3) - C(2)^2/(4*C(1));
%      
%     estimatedTau= (i+midValEstParbolicFit-1)*sSimParams.sResolution.timeBinSuperResolution-sSimParams.sResolution.timeBin/2;
%     
%     val = estimatedTau %value to find
%     tmp1 = abs(sTargetsParams.tau-val);
%     [val2 idx] = min(tmp1); %index of closest value

    %% check for points not on the grid
%     g=(-0.5:0.05:0.5);
%     estimatedTau= (i+g-1)*sSimParams.sResolution.timeBinSuperResolution - sSimParams.sResolution.timeBin/2;
%     estimatedTauVec   = exp(-1j * 2 * pi / sSimParams.T * kVec * estimatedTau);
%     estimatedSRVec= exp(1j * 2 * pi *freqShiftRatio * singleTransmitterBw*((0:(sSimParams.M-1)).'- ((sSimParams.M-1)/2))*estimatedTau);%
%     estimatedSRVec=repmat(estimatedSRVec,sSimParams.N,1);
%     
%     tmpSr=[];
%     for k=1:length(estimatedTau)
%         tmpSr = [tmpSr abs(estimatedTauVec(:,k)'*(R*diag(B_tag(:,j))*estimatedSRVec(:,k)))];
%         
%     end
% 
%     
%     
%     [val2 idx] = max(tmpSr);
%     ii=i+g(idx);
%     tmpSr(idx);
%     
%     estimatedTau2= (ii-1)*sSimParams.sResolution.timeBinSuperResolution-sSimParams.sResolution.timeBin/2
%     
%     estimatedTau= (ii-1)*sSimParams.sResolution.timeBinSuperResolution-sSimParams.sResolution.timeBin/2;
%     estimatedTauVec          = exp(-1j * 2 * pi / sSimParams.T * kVec * estimatedTau);
%     estimatedSRVec= exp(-1j * 2 * pi *freqShiftRatio * singleTransmitterBw*((0:(sSimParams.M-1))- ((sSimParams.M-1)/2))*estimatedTau);%
%     estimatedSRVec=repmat(estimatedSRVec,1,sSimParams.N);
%     aa=sSimParams.sResolution.timeBinSuperResolution;
%     i;
%     es=[es estimatedTau2];
% %     
    


    %%

    
    
    
    
    %estimatedTheta=-1+2*(j-1)/M;
    
    ii=i;
    
    if sSimParams.fixAppertreDelays
        estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:))*diag(B2(j,:));
    else
        estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:));
    end
    %toc
    
    Supp(t,3) = i;
    Supp(t,4) = j;
    
    
    SuppMat=[SuppMat vec(estimatedMatrix)];
    YVec=vec(Y);
    
    
    %x_t = YVec\SuppMat
    x_t = SuppMat\YVec;
    
    absx_t=abs(x_t);
    
    %% a
%     if t>=2
%         changeIndication= abs(pre_Xt-x_t(1:end-1));%./abs(pre_Xt);
%         [val, ind] = max(changeIndication)
%         if val>0.05
%             Supp(ind,:)=[] ; %get rid of element
%             x_t(ind)=[];
%             SuppMat(:,ind)=[];
%             sumGetRid=sumGetRid+1
%             t=t-1;
%         end 
%     end
%     pre_Xt=x_t;
% 
%     if t==5 && sumGetRid==0 ;
%         ind=[4 5];
%         
%             Supp(ind,:)=[] ; %get rid of element
%             x_t(ind)=[];
%             SuppMat(:,ind)=[];
%             sumGetRid=sumGetRid+1
%             t=t-2;
%         
%     end
%     pre_Xt=x_t;
    
    
    
    
    %% g
    
%     D_t = zeros(t,t);
%     d_t = zeros(t,1);
%     
%     m = 1;
%     while m<=t
%         r = 1;
%         while r<=t
%             D_t(m,r) = B(Supp(r,2),:)*B_tag(:,Supp(m,2))*A_tag(Supp(m,1),:)*A(:,Supp(r,1));
%             r = r+1;
%         end
%         m = m+1;
%     end
%     
%     
%     p = 1;
%     while p<=t
%         d_t(p,1) = B(Supp(p,2),:)*Y'*A(:,Supp(p,1));
%         p = p+1;
%     end
%     
%     x_t = D_t\d_t;
    
    
    % 5. Compute new residual
%     % -----------------------
%     tmp2 = 0;
%     m = 1;
%     while m<=t
%         tmp2 = tmp2 + x_t(m)*A(:,Supp(m,1))*B(Supp(m,2),:);
%         m = m+1;
%     end
    
    
    

%     if t==2
%         R = Y-2*reshape(tmp2,size(Y));
%     end
        
    
    % 6. Increment t and return to step 2
    % -----------------------------------
    
    
%     if t==numIters
%         ind=find( abs(x_t)>0.1);
%         SuppMat=SuppMat(:,ind);
%         x_t=x_t(ind);
%         t=length(ind); 
%     end
    
%         if t>numIters
%             ind=find( abs(x_t)>0.1);
%             SuppMat=SuppMat(:,ind);
%             x_t=x_t(ind);
%             Supp=Supp(ind,:);
%             t=length(ind); 
%         end

    tmp2=SuppMat*x_t;
    R = Y-reshape(tmp2,size(Y));
    norm(vec(R));
    t = t+1;
    

end
% ee1=sort(sTargetsParams.tau);
% abs(diff(sort(sTargetsParams.tau)));
% ee2=sort(es');
% %a=[ee1 ee2]
% 
X = zeros(N,M);
X_SR = cell(N,M);
% Supp(:,3)-Supp(:,1);
p = 1;


while p<=numIters
    X(Supp(p,1), Supp(p,2)) = x_t(p);
    X_SR{Supp(p,1),Supp(p,2)}=[Supp(p,3) Supp(p,4)];
    p = p+1;
end
%X   = conj(X); % < TODO: check why we use hermitic A,B, and not just TRANSPOSED
%x_t = conj(x_t);