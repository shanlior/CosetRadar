function [targets_david] = coset_nyquist(g, x,targets) % SolveAlgorithm=1=FISTA; 2=OMP

N = g.CS.N_t;
P = g.P;
Q = g.Q;
Pchopped = P - Q + 1;

Ci = g.Ci;

%k = g.CS.kappa;
%m_p = g.m_p;

% Xampling
% x, X, C: Dimensions: 1=Channel, 2=Time, 3=Bucket
H_kappa = get_H_empiric(g);
C = zeros(length(Ci), length(g.CS.kappa), P - Q + 1); % output fourier matrix (measurements)
X = fft(x,[],2);
kappa_indexes = mod(g.CS.kappa, size(X,2));

for i=1:size(X,1)
    for b=1:size(X,3)
        C(i,:,b) = X(i,kappa_indexes+1,b);
        if g.CS.normalize_H_with_division
%             C(i,:,b) = C(i,:,b) ./ H_kappa.'; %% CHECK: transpose is needed?
            C(i,:,b) = C(i,:,b); 
        else
            C(i,:,b) = C(i,:,b) .* H_kappa.';
        end
    end
end

%% self check of our equation 
% 
% %size(g.h)
% H_spectra=fft(g.h,size(x,2));
% %size(H_spectra)
% %H_spectra = g.H_spectra; 
% %size(H_spectra) 
% 
% %size(x,2)
% 
% x_check = zeros(size(x,2),size(x,3));
% 
% for  k=1:size(x_check,1)
%     for b=Q:P
%         for l=1:g.L
%             x_check(k,b-Q+1) = x_check(k,b-Q+1) + H_spectra(k)*targets.a(l) * exp(1j*2*pi*g.Ci(1)*floor(targets.t(l)/g.tau)/P) * ...
%                 exp(-2*1j*pi*(k-1)*mod(targets.t(l),g.tau)/g.tau)*exp(-1j*2*pi*(b-1)*(g.tau*targets.f(l)+g.Ci(1)/P));
%         end
%     end
% end
% DebugArray = squeeze(X(1,:,:));
% debug.a = DebugArray;
% debug.check = x_check;
% % max(max(abs(x_check-squeeze(X(1,:,:)))))
% max(max(abs(x_check-squeeze(X(1,:,:)))));


%% Focusing Cell
%clear X;

% Creates the Vandermonde matrices

% Let's annotate the matrices Y=AXB

% A matrix creation - Channels(Ci) x Q(Ambiguity factor)

% A = exp(2j*pi*Ci'*(0:Q-1)/P);
% % B matrix creation - size N x N
% % represents the delays
% 
% B = exp(-2j*pi*(0:N-1)'*(0:N-1)/N);
% % Tau_l = n * (tau/N)
% 
% targets.a = zeros(g.L,1);
% targets.t = zeros(g.L,1);
% targets.f = zeros(g.L,1);
% targets.q = zeros(g.L,1);
% disp('Running Focusing Procedure');
% shiftMatrix = zeros(size(C));
% for i = 1:size(C,1)
%     for b = 1:size(C,3)
%         shiftMatrix(i,:,b) = exp(2j*pi*b*Ci(i)/P);
%     end
% end
% % Shift by Ci
% % CHECK: if it should be done inside the loop
% shiftedC = C .* shiftMatrix;
% 
% for l=1:g.L
%     disp(strcat('Running Focusing Iteration: ', num2str(l)));
%     
%     Psi = fft(shiftedC, g.CS.N_f, 3); % Can be [] instead of g.CS.N_f
%     %Psi = fft(C, g.CS.N_f, 3); % Doppler focusing
%     % N_f is 200
%     
%     for f=1:g.CS.N_f
%         disp(strcat('Running OMP Iteration: ', num2str(f),' out of ',num2str(g.CS.N_f)));
%        [map, Supp, x_t] = OMPmatrix(Psi(:,:,f),A,B,1);
% %         map = FISTA(Psi(:,:,f),A,B,Q,N);
% 
%         map_power = abs(map).^2;
%         % debug
%         [maxRows, bestRow] = max(abs(map));
%         [~, bestCol] = max(maxRows);
%         bestRow = bestRow(bestCol);
%         totMAX = max(max(abs(map)));
%         % end of debug
% %         local_max = ordfilt2(map_power, 9, ones(3));
% 
%         %assert(length(Supp) == g.L);
%         %assert(all(abs(time_corrections) <= 0.5 + 1e-4));
%         if max(max((abs(map)))) > abs(targets.a(l))
%             [maxRows, bestRow] = max(abs(map));
%             [~, bestCol] = max(maxRows);
%             bestRow = bestRow(bestCol);
%             t_index = bestCol;
%             q_index = bestRow;
%             t_index_fix = 0;
% %             if local_max(q_index , t_index) == map_power(q_index, t_index)
% %                         checkTau = 1;
% %                 for i=1:l
% %                     for deltaTau=-(Q-1):(Q-1)
% %                         if (t_index + deltaTau * round(g.tau * g.Fs)  == targets_idx(l))
% %                             checkTau = 0;
% %                         end
% %                     end
% %                 end
% %                 if (~checkTau)
% %                     continue;
% %                 end
% %             end
%             if t_index > 1 && t_index < size(map_power, 2)
%                 t_power = map_power(q_index, t_index+(-1:1)).';
%                 Q1 = [.5 -1 .5 ; -.5 0 .5 ; 0 1 0]; % for parabola interpolation
%                 a = Q1 * t_power;
%                 assert(a(1) < 0);
%                 t_index_fix = -a(2) / 2 / a(1);
%             end
%        
%             
%             targets.a(l) = map(bestRow, bestCol);
%             targets.q(l) = bestRow - 1;
% %              targets.t(l) = ((t_index - 1 + t_index_fix) * sample_SubNyquist_factor + time_bias_fix) * g.tau/N;
%             time_bias_fix = 0.4;
%             targets.t(l) = (t_index - 1 + t_index_fix + time_bias_fix) * g.CS.delta_t;
%             %targets.t(l) = (bestCol - 1) * g.CS.delta_t;
%             targets.f(l) = (f-1) / g.CS.N_f / g.tau;
%         end
%     end
%     C_temp = zeros(size(shiftedC));
%     % subtract estimated target response
%     targets.a(l) = targets.a(l)/ (P - Q + 1);
% %     for i = 1:size(shiftedC,1)
% %         C_temp(i,:,:) = targets.a(l) * exp(-1j*2*pi*g.CS.kappa*targets.t(l)/g.tau) * ...
% %             exp(1j*2*pi*targets.f(l)*(g.Q:g.P)*g.tau) * exp(1j*2*pi*Ci(i)*targets.q(l)/g.P); %/g.tau; %CHECK: divide by Tau 
% %         % CHECK: exp(-1j*2*pi*targets.f(l)*(g.Q:g.P)*g.tau)
% %     end
%     x_cut = zeros(size(A,2),size(B,1));
%     x_cut(q_index,t_index)=targets.a(l);
% 
%     for i = 1:size(shiftedC,3)
%         C_temp(:,:,i) = A *x_cut * exp(-j*2*pi*i*ones(N,1)*(0:N-1)/N);
%         % CHECK: exp(-1j*2*pi*targets.f(l)*(g.Q:g.P)*g.tau)
%     end
%     %disp(max(max(max(abs(C_temp)))))
%     %disp(max(max(max(abs(shiftedC)))))
%     shiftedC = shiftedC - C_temp.*shiftMatrix;
%     
% end
% % disp(targets.q)
% % disp(targets.t)
% % disp(targets.f)
% targets.real_t = targets.t + g.tau * targets.q;



%% fdf

% d
% map_power = abs(map).^2;
% 
% 
% % Peak detection
% indices.t = [];
% indices.f = [];
% 
% sample_SubNyquist_factor = 1;
% [~,map_power_indexes] = sort(map_power(:),'descend');
% targets.a = zeros(g.L, 1);
% targets.t = zeros(g.L, 1);
% targets.f = zeros(g.L, 1);
% targets_idx = zeros(g.L, 1);
% l = 0;
% index_count = 1;
% local_max = ordfilt2(map_power, 9, ones(3));
% Q1 = [.5 -1 .5 ; -.5 0 .5 ; 0 1 0]; % for parabola interpolation
% while l < g.L && index_count <= length(map_power_indexes)
%     [t_index,f_index] = ind2sub(size(map_power), map_power_indexes(index_count));
%     if local_max(t_index, f_index) == map_power(t_index, f_index)
%                 checkTau = 1;
% %         for i=1:l
% %             for deltaTau=-(Q-1):(Q-1)
% %                 if (t_index + deltaTau * round(g.tau * g.Fs)  == targets_idx(l))
% %                     checkTau = 0;
% %                 end
% %             end
% %         end
% %         if (~checkTau)
% %             continue;
% %         end
%         l = l + 1;
%         targets_idx(l) = t_index;
%         targets.a(l) = map(t_index, f_index);
%         time_bias_fix = 0.4;
%         
%         % time interpolation
%         t_index_fix = 0;
%         if t_index > 1 && t_index < size(map_power, 1)
%             t_power = map_power(t_index+(-1:1), f_index);
%             a = Q1 * t_power;
% %             assert(a(1) < 0);
%             t_index_fix = -a(2) / 2 / a(1);
%         end
%         targets.t(l) = ((t_index - 1 + t_index_fix) * sample_SubNyquist_factor + time_bias_fix) * g.tau/N;
%         %assert(targets.t(l)>=0 && targets.t(l)<=g.tau*Q);
%         
%         % frequency interpolation
%         f_index_fix = 0;
%         if f_index > 1 && f_index < size(map_power, 2)
%             f_power = map_power(t_index, f_index+(-1:1));
%             a = Q1 * f_power.';
% %             assert(a(1) < 0);
%             f_index_fix = -a(2) / 2 / a(1);
%         end
%         targets.f(l) = (f_index - 1 + f_index_fix) /(P*g.tau);
% %         assert(targets.f(l)>=0 && targets.f(l)<=1/g.tau);
%         
%         if g.mf.debug_plot
%             h = plot(targets.t(l), targets.f(l), 'ro');
%             set(h,'MarkerSize',10);
%         end
%     end
%     index_count = index_count + 1;
%     indices.t = [indices.t ; t_index];
%     indices.f = [indices.f ; f_index];
% end
%     [indices.t' ; indices.f']';
    
%% Kronecker

% Changes: h (pulse), P size, h_BW
% x, X, C: Dimensions: 1=Channel, 2=Time, 3=Bucket


% Creates the Vandermonde matrices

% Let's annotate the matrices Y=AXB

% A matrix creation - size N x (Q*N)
% represents the delays

A = zeros(size(C,1),N,Q*N);
% A2 = zeros(size(C,1),N,Q*N);


A1 = exp(-j*2*pi*(0:N-1)'*(0:N-1)/N);
for i=1:size(A1,2)
    A1(:,i) = A1(:,i) .* H_kappa;
end
% disp STARTa
% tic
for c = 1:size(C,1)
    for q=1:Q
        A(c,:,(1+(q-1)*N):(q*N)) = A1 * exp(j*2*pi*Ci(c)*(q-1)/P);
    end
end
% disp 'start'
% for c = 1:size(C,1)
%     for n1=1:N
%         for q=1:Q
%             for n2=1:N
%                 A2(c,n1,(q-1)*N+n2) = exp(j*2*pi*Ci(c)*(q-1)/P -j *2*pi*(n1-1)'*(n2-1)/N);
%             end
%         end
%     end
% end
% disp end


% disp A
% toc
% ersds
% B matrix creation - size P x (P-Q)

B = zeros(size(C,1),P,P-Q+1);
% tic
for c = 1:size(C,1)
    for p = 1:P
        for b = Q:P
            %B(p,b-Q+1) = exp(-(j*2*pi/P)*((p-1)*(b-1) - Ci*(b-1)));
            B(c,p,b-Q+1) = exp(-2j*pi*(b-1)*((p-1)/P + Ci(c)/P));
        end
    end
end
% disp B
% toc
[X, R, X_SR,Supp] = SolveOmpDavid(C, A, B, g.L);
targets_david.t = squeeze(Supp(:,1));
targets_david.f = squeeze(Supp(:,2));

%% kron_final
% krons = [];
% Cvec = [];
% for c = 1:size(C,1)
%     krons = [krons ; kron(transpose(squeeze(B(c,:,:))),squeeze(A(c,:,:)))]; % solutions
%     Cvec = [Cvec ; reshape(squeeze(C(c,:,:)),size(C,2)*size(C,3),1)]; % measures
% end
% t_kron = zeros(g.L,1);
% f_kron = zeros(g.L,1);
% a_kron = zeros(g.L,1);
% map_rec = pinv(krons) * Cvec;
% map_rec_reshape = reshape(map_rec,Q*N, P);
% map = reshape(abs(map_rec),Q*N, P);
% for i = 1:g.L
%     [maxRows, bestRow] = max(map);
%     [~, bestCol] = max(maxRows);
%     bestRow = bestRow(bestCol);
%     while any(find(abs(bestRow - t_kron) <= 1))
%             map(bestRow,bestCol) = 0;
%             [maxRows, bestRow] = max(map);
%             [~, bestCol] = max(maxRows);
%             bestRow = bestRow(bestCol);
%     end
%     t_kron(i) = bestRow;
%     f_kron(i) = bestCol;
%     a_kron(i) = map_rec_reshape(bestRow,bestCol);
%     map(bestRow,bestCol) = 0;
%     if bestRow < size(map,1)
%         map(bestRow+1,bestCol) = 0;
%     end
%     if bestRow > 1
%         map(bestRow-1,bestCol) = 0;
%     end
% end
% targets_kron.t = t_kron;
% targets_kron.f = f_kron;
% targets_kron.a = a_kron;
% 
% % end Kron


%% OMP VEC 
% 
% % Creates the Vandermonde matrices
% 
% % Let's annotate the matrices Y=AXB
% 
% % A matrix creation - size N x (Q*N)
% % represents the delays
% 
% A = zeros(size(C,1),N,Q*N);
% tmp = exp(-j*2*pi*(0:N-1)'*(0:N-1)/N);
% for c = 1:size(C,1)
%     for q=1:Q
%         A(c,:,(1+(q-1)*N):(q*N)) = tmp * exp(j*2*pi*Ci(c)*(q-1)/P);
%     end
% end
% 
% % B matrix creation - size P x (P-Q)
% 
% B = zeros(size(C,1),P,P-Q+1);
% for c = 1:size(C,1)
%     for p = 1:P
%         for b = Q:P
%             %B(p,b-Q+1) = exp(-(j*2*pi/P)*((p-1)*(b-1) - Ci*(b-1)));
%             B(c,p,b-Q+1) = exp(-2j*pi*(b-1)*((p-1)/P + Ci(c)/P));
%         end
%     end
% end
% krons = [];
% Cvec = [];
% for c = 1:size(C,1)
%     krons = [krons ; kron(transpose(squeeze(B(c,:,:))),squeeze(A(c,:,:)))]; % solutions
%     Cvec = [Cvec ; reshape(squeeze(C(c,:,:)),size(C,2)*size(C,3),1)]; % measures
% end
% t_OMP_Vec = zeros(g.L,1);
% f_OMP_Vec = zeros(g.L,1);
% 
% [Xr, Supp, ResNorm, NormResvsSol] = RunOMP_Unnormalized(Cvec,krons,4, 0, 0.01, false);
% map = reshape(abs(Xr),Q*N, P);
% [targets_OMP_vec.t,targets_OMP_vec.f]=find(map);



end
