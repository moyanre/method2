function [L_Mat_Loni,Z_Mat_Loni] = fn_L_Z_cIBP_Trans(L_Mat_Lana,Z_Mat_Lana,alpha,beta,Q,time_step)



if time_step == 1
    
            K_i = poissrnd(alpha); 

                weit = (beta*ones(1,Q))/(Q*beta);
                first_row = zeros(1,K_i);

                for d = 1:K_i
                    first_row(d) = discretesample(weit,1);
                end

                L_Mat_Loni = first_row;

                
                
                %%%%%%%%  Z part 
                
                        len = length(L_Mat_Loni);

                        Z_Mat_Loni = zeros(1,len);
                        for f = 1:len
                            hee = L_Mat_Loni(f);
                            
                            if hee == 0
                                Z_Mat_Loni(f) = 0;
                            else
                                prob = (1/(hee+1))*ones(1,hee+1); 
                                chosee = discretesample(prob,1);
                                Z_Mat_Loni(f) = chosee - 1;
                            end
                        end
    
    
    
else
   
    non_zero_col = zeros(1,size(L_Mat_Lana,2));
    
    for e = 1:size(L_Mat_Lana,2)
        col_j = L_Mat_Lana(:,e);
        
        if sum(col_j) > 0
           non_zero_col(e) = 1; 
        end
    end
    
    K_plus = sum(non_zero_col);
    
    new_row_part_1 = zeros(1,K_plus);
    
    
    for k = 1:K_plus
        col_k = L_Mat_Lana(:,k);
        m_k = length(find(col_k ~= 0));
        
        
        prob_Q = zeros(1,Q);
        
        for q = 1:Q
            
            m_kq = length(find(col_k == q));
            
            prob_Q(q) = (m_k/time_step) * ((beta + m_kq)/((beta*Q) + m_k));
            
        end
        
        
        %m_k1 = length(find(col_k == 1));
        %m_k2 = length(find(col_k == 2));
        %m_k3 = length(find(col_k == 3));
        %p1 = (m_k/time_step) * ((betas(1) + m_k1)/(sum(betas) + m_k));
        %p2 = (m_k/time_step) * ((betas(2) + m_k2)/(sum(betas) + m_k));
        %p3 = (m_k/time_step) * ((betas(3) + m_k3)/(sum(betas) + m_k));
        %p0 = 1 - (p1 + p2 + p3);
        
        p0 = 1 - sum(prob_Q);
        prob_all = [p0 prob_Q]; 
        
        
        val = discretesample(prob_all,1);
        new_row_part_1(k) = val - 1;
        
    end
    
    
    K_i_new = poissrnd(alpha/time_step);
    
            if K_i_new == 0
                
                L_Mat_Loni = [L_Mat_Lana;new_row_part_1];

            else
                
                
                new_row_part_2 = zeros(1,K_i_new);
                      weit_2 = beta/sum(beta);
                       for v = 1:K_i_new
                           new_row_part_2(v) = discretesample(weit_2,1);
                       end
                
                
                new_row = [new_row_part_1 new_row_part_2];


                padding_zeros = zeros(size(L_Mat_Lana,1),K_i_new);

                L_Mat_Loni = [L_Mat_Lana  padding_zeros; new_row];

            end
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%% Z_part
                        new_row_L = L_Mat_Loni(end,:);
                        len = length(new_row_L);

                        Z_new_row = zeros(1,len);
                        for f = 1:len
                            hee = new_row_L(f);

                                if hee == 0
                                   Z_new_row(f) = 0;
                                else
                                    prob = (1/(hee+1))*ones(1,hee+1);
                                    chosee = discretesample(prob,1);
                                    Z_new_row(f) = chosee - 1;  
                                end


                        end


                           ro_la = size(Z_Mat_Lana,1);
                           co_la = size(Z_Mat_Lana,2);

                        if co_la == len

                             Z_Mat_Loni = [Z_Mat_Lana ; Z_new_row];

                        else

                             Z_Mat_Loni = [Z_Mat_Lana zeros(ro_la,len-co_la); Z_new_row];

                        end
 
  
end




end


            




