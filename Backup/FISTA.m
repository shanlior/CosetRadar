function X = FISTA(Y,A,B,N,M)
    beta=0.9;
    L=1000000;
    lambda_bar=1e-7;
    lambda=1e-5;
    t_k=1;
    t_k_m1=1;
    X=zeros(N,M);
    X_k_m1=X;
    
    for i=1:30
        Z=X+((t_k_m1-1)/t_k)*(X-X_k_m1);
        grad_f_z=(A'*(A*(Z)*B-Y)*B');
        U=Z-(1/L)*grad_f_z;
        disp(['Iteration #:' num2str(i)]);
        X=sign(U).*(abs(U)-lambda/L).*(abs(U)>lambda/L);
        X_k_m1=X;
        t_k_m1=t_k;
        t_k=(1+sqrt(4*t_k*t_k+1))/2;
        lambda=max(beta*lambda,lambda_bar);
    end


    