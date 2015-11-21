
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
 
    if sSimParams.fixAppertreDelays
        tmp=(A2'.*(A1'*R))*(B.*B2)';
    else
        tmp=(A2'.*(A1'*R))*B';
    end
    %toc
    
    %debug mode
   [~,i_vec] = max(tmp(:));
   [i, j] = ind2sub(size(tmp),ind);
    
    % 3. Augment index set
    % --------------------
   Supp(t,2) = j;  %sin
   Supp(t,1) = i;  %range 
   
    
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
    YVec=Y(:);
    
    
    x_t = SuppMat\YVec;
    
    tmp2=SuppMat*x_t;
    R = Y-reshape(tmp2,size(Y));
    norm(vec(R));
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