clear all; clc; close all; rng('default'); rng(100);  
load('counts_alt'); load('counts_total')            
                    
Data1 = counts_alt; N_st_Matrix1 = counts_total;
QQ = 3; Q_vec2 = QQ*ones(size(N_st_Matrix1,1),1);
 
rep = 5; Data = repmat(Data1,rep,1);
N_st_Mat = repmat(N_st_Matrix1,rep,1); Q_vec = repmat(Q_vec2,rep,1);
S = size(Data,1); T = size(Data,2);
%%%%% initialization
        % parameters of cIBP
          alpha = 0.8; beta = 2;
        % number of particles
          N = 500; 
        % Initialize p_o
          a_00 = 1; b_00 = 100; p_o_particles = betarnd(a_00,b_00,[1,N]);
        % Initialize theta_ot, t = 1,...,T
          a_0 = 0.1; b_general = 1; theta_o_particles = gamrnd(a_0,b_general,[T,N]);
          
        % Initialize phi_t t = 1,...,T 
          b_phi = 1; mean_N_st = median(N_st_Matrix1);
          a_phi = mean_N_st*0.5;
          
          phi_particles = zeros(T,N);
          for t = 1:T
              phi_particles(t,:) = gamrnd(a_phi(t),1/b_phi,[1,N]);
          end
        % parameters of theta_ct
          a = 1;
        % Initialize matrix L and Z
          L_init = 0; Z_init = 0;  
          para_frac = 0.01; %0.03          
for s = 1:S
    Data_s = Data(s,:); N_st_s = N_st_Mat(s,:); s
    if s == 1
            for n = 1:N
                L_lana_n = L_init; Z_lana_n = Z_init;
                Q = Q_vec(s);
                [L_loni_n, Z_loni_n] = fn_L_Z_cIBP_Trans(L_lana_n,Z_lana_n,alpha,beta,Q,s);
                L_loni{n} = L_loni_n; Z_loni{n} = Z_loni_n;
                C_here_n = size(L_loni_n,2);
                no_theta_n = C_here_n*T;
                theta_ct_n = gamrnd(a,b_general,[no_theta_n,1]);
                %%% rejuvenate theta_ct 
                theta_ct_n = normrnd(theta_ct_n,para_frac*theta_ct_n);
                theta_ct_ALL{n} = theta_ct_n;  %%% all theta_ct
                %%%% rejuvenate p_o
                p_o_particles(n) = normrnd(p_o_particles(n),para_frac*p_o_particles(n));
                po_n = p_o_particles(n);
                %%%% rejuvenate theta_o
                theta_o_particles(:,n) = normrnd(theta_o_particles(:,n),para_frac*theta_o_particles(:,n));
                theta_o_n_all_T = theta_o_particles(:,n);
                %%%% rejuvenate phi's
                phi_particles(:,n) = normrnd(phi_particles(:,n),para_frac*phi_particles(:,n));
                phi_n = phi_particles(:,n);
                %%% Calculate the weight for particle n
                prob_n_nn = vpa(zeros(1,T)); prob_n_N = vpa(zeros(1,T));
                       for t = 1:T
                           n_st = Data_s(t); NN_st = N_st_s(t); phi_t = phi_n(t);
                           theta_o_t = theta_o_n_all_T(t);
                         %%% take out the theta_ct
                           bee = (t-1)*C_here_n + 1; eed = t*C_here_n; 
                           theta_ct = theta_ct_n(bee:eed);
                           all_theta = [theta_o_t; theta_ct]; 
                           w_ct = all_theta/sum(all_theta); 
                           lo_and_L = [2 L_loni_n(s,:)]; 
                           po_and_z = [2*po_n Z_loni_n(s,:)]; 
                           M_st = lo_and_L*w_ct;
                           p_st = (po_and_z*w_ct)/M_st;
                           prob_n_nn(t) = my_binopdf(n_st,NN_st,p_st);
                           prob_n_N(t) = my_poisspdf(NN_st,(0.5*phi_t*M_st));
                       end
                          inter1 = exp(sum(log([prob_n_nn prob_n_N])));
                          PROB(n) = inter1;  
            end
            %%% normalize weights
                wt = PROB/sum(PROB); wt = double(wt);
            %%% resample  %%% rr = randsample([1:N],N,true,wt), new = old{rr}
                 rr = randsample([1:N],N,true,wt);
                 p_o_particles = p_o_particles(rr);
                 theta_o_particles = theta_o_particles(:,rr);
                 phi_particles = phi_particles(:,rr);
                 L_loni = L_loni(rr); Z_loni = Z_loni(rr);
                 theta_ct_ALL = theta_ct_ALL(rr);
    else 
        Z_lana = Z_loni; L_lana = L_loni; 
            for n = 1:N
                L_lana_n = L_lana{n}; Z_lana_n = Z_lana{n}; Q = Q_vec(s);
                [L_loni_n, Z_loni_n] = fn_L_Z_cIBP_Trans(L_lana_n,Z_lana_n,alpha,beta,Q,s);
                L_loni{n} = L_loni_n; Z_loni{n} = Z_loni_n;
                C_lana = size(L_lana_n,2); C_loni = size(L_loni_n,2); diff = C_loni - C_lana;
                %%% rejuvenate p_o
                p_o_particles(n) = normrnd(p_o_particles(n),para_frac*p_o_particles(n));
                po_n = p_o_particles(n);
                %%% rejuvenate theta_o
                theta_o_particles(:,n) = normrnd(theta_o_particles(:,n),para_frac*theta_o_particles(:,n));
                theta_o_n_all_T = theta_o_particles(:,n);
                %%% rejuvenate theta_ct
                theta_ct_ALL{n} = normrnd(theta_ct_ALL{n},para_frac*theta_ct_ALL{n});
                theta_ct_n = theta_ct_ALL{n};  %%% all theta_ct
                %%%% rejuvenate phi's
                phi_particles(:,n) = normrnd(phi_particles(:,n),para_frac*phi_particles(:,n));
                phi_n = phi_particles(:,n);
                  if diff == 0  %%% new column(s) of Z is/are not created.
                            prob_n_nn = vpa(zeros(1,T)); prob_n_N = vpa(zeros(1,T));
                               for t = 1:T 
                                   n_st = Data_s(t); NN_st = N_st_s(t);
                                   phi_t = phi_n(t); theta_o_t = theta_o_n_all_T(t);
                                 %%% take out the theta_ct
                                   bee = (t-1)*C_lana + 1; eed = t*C_lana; 
                                   theta_ct = theta_ct_n(bee:eed);
                                   all_theta = [theta_o_t; theta_ct]; 
                                   w_ct = all_theta/sum(all_theta); 
                                   lo_and_L = [2 L_loni_n(s,:)]; 
                                   po_and_z = [2*po_n Z_loni_n(s,:)]; 
                                   M_st = lo_and_L*w_ct;
                                   p_st = (po_and_z*w_ct)/M_st;
                                   prob_n_nn(t) = my_binopdf(n_st,NN_st,p_st);
                                   prob_n_N(t) = my_poisspdf(NN_st,(0.5*phi_t*M_st));
                               end      
                  else
                                 prob_n_nn = vpa(zeros(1,T)); prob_n_N = vpa(zeros(1,T));
                                 theta_ct_new_matrix = zeros(C_loni,T);
                               for t = 1:T
                                   n_st = Data_s(t); NN_st = N_st_s(t); phi_t = phi_n(t);
                                   theta_o_t = theta_o_n_all_T(t);
                                 %%% take out the theta_ct
                                   bee = (t-1)*C_lana + 1; eed = t*C_lana; 
                                   theta_ct = theta_ct_n(bee:eed);
                                   theta_ct_loni = [theta_ct_n(bee:eed); gamrnd(a,b_general,[diff,1])];
                                   theta_ct_new_matrix(:,t) = theta_ct_loni;
                                   all_theta = [theta_o_t; theta_ct_loni]; 
                                   w_ct = all_theta/sum(all_theta);
                                   lo_and_L = [2 L_loni_n(s,:)];
                                   po_and_z = [2*po_n Z_loni_n(s,:)]; 
                                   M_st = lo_and_L*w_ct; p_st = (po_and_z*w_ct)/M_st;
                                   prob_n_nn(t) = my_binopdf(n_st,NN_st,p_st);
                                   prob_n_N(t) = my_poisspdf(NN_st,(0.5*phi_t*M_st));
                               end
                                 theta_ct_ALL{n} = theta_ct_new_matrix(:);
                  end                              
                                        inter3 = sum(log([prob_n_nn prob_n_N]));
                                        inter4 = exp(vpa(inter3)); PROB(n) = inter4;    
            end
            %%% normalize weights
                wt = PROB/sum(PROB); wt = double(wt);
                rr = randsample([1:N],N,true,wt); p_o_particles = p_o_particles(rr);
                theta_o_particles = theta_o_particles(:,rr);
                phi_particles = phi_particles(:,rr);
                theta_ct_ALL = theta_ct_ALL(rr);
                L_loni = L_loni(rr); Z_loni = Z_loni(rr); 
    end    
end
wt = (1/N)*ones(1,N); S_orig = S/rep;
for n = 1:N
L_loni_now{n} = L_loni{n}(S - S_orig + 1:S,:);
Z_loni_now{n} = Z_loni{n}(S - S_orig + 1:S,:);
end
%%%%%%%  POINT ESTIMATES %%%%%%%%%%
        ALL_C = zeros(1,N);
        for n = 1:N
            ALL_C(n) = size(L_loni_now{n},2); 
        end   
       unique_C = unique(ALL_C);
       prob_C = zeros(1,length(unique_C));
       len_Cs = zeros(1,length(unique_C));
       for c = 1:length(unique_C)
           prr = unique_C(c); ind_c = find(ALL_C == prr);
           prob_C(c) = length(ind_c)/N; len_Cs(c) = length(ind_c);
       end
       [vaa, indd] = max(len_Cs); C_est = unique_C(indd)
       %%%% Plot of C
uppe = max(unique_C) + 1; lowe = min(unique_C) - 1;
x_min = 0; x_max = uppe + 5;
pre_x = 0:lowe; post_x = uppe:x_max;
xx = [pre_x unique_C post_x]; 
yy = [zeros(1,lowe+1) prob_C zeros(1,x_max-max(unique_C))]; 
figure;
plot(xx,yy)
xlabel('C')
ylabel('P(C|Y)')
ylim([0 1])
h = vline(C_est,'g','The Answer');
           
     %%%%% pick out all particles (L,Z,W,phi,p) that belong to this C-star
     ind_C_est = find(ALL_C == C_est);
    
    picked_L_loni = L_loni_now(ind_C_est);
    picked_Z_loni = Z_loni_now(ind_C_est);
    picked_p_o_particles = p_o_particles(ind_C_est); 
    picked_theta_o_particles = theta_o_particles(:,ind_C_est);
    picked_theta_ct_ALL = theta_ct_ALL(ind_C_est);
    picked_phi_particles = phi_particles(:,ind_C_est);
    picked_wt = wt(ind_C_est); 
     %%%%% determine point estimate of L   i.e L_est
   NNN = length(ind_C_est); arg_vector = zeros(1,NNN);
   all_permutations = perms(1:C_est);
     for g = 1:NNN
         L_prime = picked_L_loni{g}; dist_L_l_prime = zeros(1,NNN);
         for l = 1:NNN
             L_l = picked_L_loni{l};
                   d_min_vec = zeros(1,size(all_permutations,1)); 
                for permu = 1:size(all_permutations,1)
                    perm_here = all_permutations(permu,:);
                    permuted_L_prime = L_prime(:,perm_here);
                    d_min_vec(permu) = sum(sum(abs(L_l - permuted_L_prime)));
                end
             dist_L_l_prime(l) = min(d_min_vec);
         end
  arg_vector(g) = sum(picked_wt.*dist_L_l_prime);
     end
     
     [val,indx] = min(arg_vector);
     %%%% point estimate of L
     L_est = picked_L_loni{indx}      
     %%%% point estimate of Z  
     Z_est = picked_Z_loni{indx}     
     %%%% point estimate of p_0
     p_o_est = picked_p_o_particles(indx);
     
     %%%% point estimate of phis
     phi_est = picked_phi_particles(:,indx);
     
     %%%% point estimate of W
     theta_o_est = picked_theta_o_particles(:,indx);
     theta_ct_ALL_est = picked_theta_ct_ALL{indx};
     olu = vec2mat(theta_ct_ALL_est,C_est)';
     ade = [theta_o_est'; olu]; sum_ade = sum(ade);
     W_est = ade./repmat(sum_ade,C_est+1,1)
   