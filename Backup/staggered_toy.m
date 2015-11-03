clear all;
close all;
clc;

N=10; % residual delay grid size and number of samples per pulse
P=10; % Doppler grid size and number of pulses
Q=2; % max number of PRIs for delay

c1=2; c2=8; c3=9;

% delay Doppler matrix
A=zeros(Q*N, P);
A(3, 5)=2;
%A(8,8)=3;
A(18, 8)=3;

% Fourier delay matrix
F_t1 =[dftmtx(N) exp(-1i*2*pi*c1/P)*dftmtx(N)];
F_t2 =[dftmtx(N) exp(-1i*2*pi*c2/P)*dftmtx(N)];
F_t3 =[dftmtx(N) exp(-1i*2*pi*c3/P)*dftmtx(N)];

% Fourier Dopler matrix
F_d=dftmtx(P);

v1=zeros(1,P);
v2=zeros(1,P);
v3=zeros(1,P);
for b=1:P
    v1(b) = exp(-1i*2*pi*(b-1)*c1/P);
    v2(b) = exp(-1i*2*pi*(b-1)*c2/P);
    v3(b) = exp(-1i*2*pi*(b-1)*c3/P);
end    
V1=repmat(v1, P, 1);
F_d1=F_d.*V1;
V2=repmat(v2, P, 1);
F_d2=F_d.*V2;
V3=repmat(v3, P, 1);
F_d3=F_d.*V3;

F_d1=F_d1(:,1:P-1);
F_d2=F_d2(:,1:P-1);
F_d3=F_d3(:,1:P-1);

X1=F_t1*A*F_d1;
x1 = reshape(X1, N*(P-1), 1);
G1=kron(transpose(F_d1), F_t1);


X2=F_t2*A*F_d2;
x2 = reshape(X2, N*(P-1), 1);
G2=kron(transpose(F_d2), F_t2);

X3=F_t3*A*F_d3;
x3 = reshape(X3, N*(P-1), 1);
G3=kron(transpose(F_d3), F_t3);

G= [G1 ; G2 ; G3];
x= [x1; x2; x3];
a=pinv(G)*x;
A_rec=reshape(a, Q*N, P);
abs(A_rec);
size(F_d1)
size(F_t1)


% b=pinv(G1)*x1;
% B_rec=reshape(b, Q*N, P);
% abs(B_rec)
% 
% Gd= [G1 ; G2];
% xd= [x1; x2];
% h=pinv(Gd)*xd;
% H_rec=reshape(h, Q*N, P);
% abs(H_rec)
