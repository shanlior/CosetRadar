
function [X, R, X_SR, Supp] = SolveOmpDavid(Y, A, B, numIters)
% solves Y = A*X*B.'
% Input:
%   Y, A, B     - matrices of equation
%   numIters    - number of algorithm iterations
% Output:
%   sResultsOmp - results by OMP algorithm


% 1. Initialization
% ------------------
% x, X, C: Dimensions: 1=Channel, 2=Time, 3=Bucket
[T,N,QN] = size(A);
[T,P,P_Q] = size(B);
% Y = [];
phi = zeros(T*P_Q,QN);
A_2d = zeros(T*N,QN);
B_2d = zeros(T*P_Q,P);
estimatedMatrix = zeros(T,N,P_Q);
for (m = 1:T)
    B_2d(((m-1)*P_Q+1):(m*P_Q),:) = squeeze(B(m, :,:)).';
end



B_2d = conj(B_2d);
R = Y; % < Residual
Supp = [];  % Support
SuppMat=[];
t = 1;
M=T;



while t<=numIters %| norm(vec(R))> 80
%     str = ['Iteration #', num2str(t)];
%     disp(str);
    
    % 2. Find the two indicies of the support
    % ---------------------------------------
    for (m = 1:T)
        phi(((m-1)*P_Q+1):(m*P_Q),:) = (squeeze(A(m,:,:))' * squeeze(R(m,:,:))).';
    end

    tmp = phi.' * B_2d;
    
   
   [~,ind] = max(tmp(:));
   [i, j] = ind2sub(size(tmp),ind);
    
    % 3. Augment index set
    % --------------------
   Supp(t,2) = j;  %sin
   Supp(t,1) = i;  %range 
   
   
    ii=i;
    
     for (m= 1:T)
         A_2d(((m-1)*N+1):(m*N),:) = squeeze(A(m,:,:));
     end
     
     for (m = 1:T)
         estimatedMatrix(m,:,:) = squeeze(A(m,:,i)).' * squeeze(B(m,j,:)).';
     end

    Supp(t,3) = i;
    Supp(t,4) = j;
    [y1,y2,y3] = size(estimatedMatrix);
    SuppMat=[SuppMat reshape(estimatedMatrix,y1*y2*y3,1)];
    [y1,y2,y3] = size(Y);
    YVec=reshape(Y,y1*y2*y3,1);

    x_t = SuppMat\YVec;
        
 
    tmp2=SuppMat*x_t;
    R = Y-reshape(tmp2,size(Y));
    
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
