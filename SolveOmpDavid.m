
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
B_2d = [];
for (iii = 1:T)
    B_2d = [B_2d ; squeeze(B(iii, :,:)).'];
%     Y = [Y ; squeeze(C(iii,:,:))];
end


B_2d = conj(B_2d);
R = Y; % < Residual
Supp = [];  % Support
SuppMat=[];
t = 1;
M=T;
sumGetRid=0;



while t<=numIters %| norm(vec(R))> 80
    str = ['Iteration #', num2str(t)];
    disp(str);
    
    % 2. Find the two indicies of the support
    % ---------------------------------------
 
    phi = [];
    for (m = 1:T)
        phi = [phi ; (squeeze(A(m,:,:))' * squeeze(R(m,:,:))).'];
    end
    phi = phi.';
    tmp = phi * B_2d;
    
%     if sSimParams.fixAppertreDelays
%         tmp=(A2'.*(A1'*R))*(B.*B2)';
%     else
%         tmp=(A2'.*(A1'*R))*B';
%     end
 
   [~,i_vec] = max(tmp(:));
   [i, j] = ind2sub(size(tmp),ind);
    
    % 3. Augment index set
    % --------------------
   Supp(t,2) = j;  %sin
   Supp(t,1) = i;  %range 

   
   
    ii=i;
    
%     if sSimParams.fixAppertreDelays
%         estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:))*diag(B2(j,:));
%     else
%         estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:));
%     end
%     %toc
%     A1 = exp(-1j * 2 * pi / sSimParams.T * kVec * tauVec);
%      A2 = exp(-1j * 2 * pi *fVec  * tauVec);
% 
%     A1 = kvec * tauVec
%     A2 = fvec *tauVec

     A_2d = [];
     for (m= 1:T)
         A_2d = [A_2d ; squeeze(A(m,:,:))];
     end
     
     
     for (m = 1:T)
         estimatedMatrix(m,:,:) = squeeze(A(m,:,i)).' * squeeze(B(m,j,:)).';
     end
%       size(A)
%      size(B)
%      size(A_2d)
% 
%      size(conj(B_2d))
%      size(tmp)

%      A(m,:,i) * conj(B_2d(m,j,:));
%      estimatedMatrix = * conj(B_2d(j,:));
%      estimatedMatrix = zeros(
%      estimatedMatrix=A1(:,i)*A2(:,i).'*diag(B(j,:))*diag(B2(j,:));

    Supp(t,3) = i;
    Supp(t,4) = j;
    [y1,y2,y3] = size(estimatedMatrix);
    SuppMat=[SuppMat reshape(estimatedMatrix,y1*y2*y3,1)];
    [y1,y2,y3] = size(Y);
    YVec=reshape(Y,y1*y2*y3,1);

    %x_t = YVec\SuppMat
    x_t = SuppMat\YVec;
    
    absx_t=abs(x_t);
    
 
    tmp2=SuppMat*x_t;
    R = Y-reshape(tmp2,size(Y));
    norm(reshape(R,y1*y2*y3,1));
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
