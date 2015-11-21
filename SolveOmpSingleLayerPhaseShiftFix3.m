
function [X, R, X_SR] = SolveOmpSingleLayerPhaseShiftFix3(Y, A1, B,A2 ,B2, numIters,sSimParams)
% solves Y = A*X*B.'
% Input:
%   Y, A1, B,A2 ,B2     - matrices of equation
%   numIters    - number of algorithm iterations
% Output:
%   sResultsOmp - results by OMP algorithm


% 1. Initialization
% ------------------
R = Y; % < Residual
Supp = [];  % Support
SuppMat=[];
t = 1;
[~,N] = size(A1);
[M,~] = size(B);

while t<=numIters %| norm(vec(R))> 80
    %     str = ['Iteration #', num2str(t)];
    %     disp(str);
    
    % 1. Find the two indicies of the support
    % ---------------------------------------
    %tic
    if sSimParams.fixAppertreDelays
        tmp=(A2'.*(A1'*R))*(B.*B2)';
    else
        tmp=(A2'.*(A1'*R))*B';
    end
    [~ , ind] =max(tmp(:));
    [i , j] = ind2sub(size(tmp) ,ind)
    %toc
    
    %debug mode
    %stem3(abs(tmp))
    %view(180,0)
    %    if t==1
    %     figure;
    %     subplot 311
    %     stem3(abs(tmp))
    %     subplot 312
    %     stem3(abs(tmp))
    %     view(90,0)
    %     subplot 313
    %     stem3(abs(tmp))
    %     view(180,0)
    %     ade=1;
    %    end
    

    % 2. Update Support
    % --------------------
    Supp(t,1) = i;  %range
    Supp(t,2) = j;  %sin
%     estimatedTau= (i-1)*sSimParams.sResolution.timeBinSuperResolution;
%     estimatedTheta=-1+2*(j-1)/M;   
    
    
    
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

    %%
    
    % 3. Update coffecients
    % -------------------------------
    %calculate matrix for the new entry
 
    if sSimParams.fixAppertreDelays
        estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:))*diag(B2(j,:));
    else
        estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:));
    end

    Supp(t,3) = i;
    Supp(t,4) = j;
    
    SuppMat=[SuppMat estimatedMatrix(:)];
    YVec=Y(:);
    
    x_t = SuppMat\YVec;
    
    %absx_t=abs(x_t);
    
    
    % 4. Compute new residual
    %     % -----------------------
    
    tmp2=SuppMat*x_t;
    R = Y-reshape(tmp2,size(Y));

    %norm(R(:));
    t = t+1;
    
end

X = zeros(N,M);
X_SR = cell(N,M);
p = 1;

while p<=numIters
    X(Supp(p,1), Supp(p,2)) = x_t(p);
    X_SR{Supp(p,1),Supp(p,2)}=[Supp(p,3) Supp(p,4)];
    p = p+1;
end
