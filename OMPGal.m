function [X, Supp, x_t] = OMPGal(Y, A, B, NumIters)

% 1. Initialization
% ------------------
R = Y;      % Residual
Supp = [];  % Support
t = 1;
A_tag = A';
B_tag = B';
[g_N,N] = size(A);
Q=N/g_N;
[M,~] = size(B);

while t<=NumIters  
    str = ['Iteration #', num2str(t)];
    disp(str);
    
    % 2. Find the two indicies of the support
    % ---------------------------------------
    
    tmp1 = abs(A_tag*R*B_tag);
  while(1) 
      flag=0;
        [~,i_vec] = max(tmp1);
        [max_val, j] = max(max(tmp1));
        i = i_vec(j);

        if (size(Supp,1)~=0)
            for rr=1:size(Supp,1)
                for kk=1:Q
                    if (Supp(rr,1)==mod(i+kk*g_N,N) || Supp(rr,1)==mod(i+kk*g_N,N)+1 || Supp(rr,1)==mod(i+kk*g_N,N)-1)
                        flag=1; 
                    end
                end
            end
        end
        if flag
             tmp1(i,j)=0;
            continue
        else
            Supp(t,2) = j;
            Supp(t,1) = i;
            break
        end
  end
    
%     
%     
%     % 3. Augment index set
%     % --------------------
%     Supp(t,2) = j;
%     Supp(t,1) = i;
%     
%     Supp_uniq = unique(Supp, 'rows');
%     for ii=1:Q-1
%         supp_tmp_i=i+ii*g.N;
%         
%      end
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
    D_t = zeros(t,t);
    d_t = zeros(t,1);
    
    m = 1;
    while m<=t
        r = 1;
        while r<=t
            D_t(m,r) = B(Supp(r,2),:)*B_tag(:,Supp(m,2))*A_tag(Supp(m,1),:)*A(:,Supp(r,1));
            r = r+1;
        end
        m = m+1;
    end
    
    
    p = 1;
    while p<=t
        d_t(p,1) = B(Supp(p,2),:)*Y'*A(:,Supp(p,1));
        p = p+1;
    end
    
    x_t = D_t\d_t;
    
    
    % 5. Compute new residual
    % -----------------------
    tmp2 = 0;
    m = 1;
    while m<=t
        tmp2 = tmp2 + x_t(m)*A(:,Supp(m,1))*B(Supp(m,2),:);
        m = m+1;
    end
    
    R = Y-tmp2;
    
    % 6. Increment t and return to step 2
    % -----------------------------------
    t = t+1;
end

X = zeros(N,M);

p = 1;
while p<=NumIters
    X(Supp(p,1),Supp(p,2)) = x_t(p);
    p = p+1;
end