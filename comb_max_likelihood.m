function [max_theta_site,max_p_array,max_alpha,max_likelihood,g_mix_sel_fold,g_mix_neu_fold] = comb_max_likelihood(selected_sfs,neutral_sfs,file_name,initial_gamma_array,boundary_array,initial_p_array,initial_theta_site,initial_lethal_array,free_alpha,neutral_alpha,zero_class,N_pop)
%selected_sfs folded, neutral_sfs folded, assumes they contain 0-class

    function [P_samp] = sample_prob(Npop,Nsamp)
        x = (0:(2*Npop-1))/(2*Npop);
        P_samp = zeros((Nsamp+1),length(x));
        for k = 0:Nsamp
        	P_samp((k+1),:) = binopdf(k,Nsamp,x);
        end
    end

    N__samp = 2*(length(selected_sfs)-1);
    P__samp = sample_prob(N_pop,N__samp);

    function [G,m_samp] = theor_SFS(gamma_array,mu_persite,Num_sites,N_pop,N_samp,P_samp)
        %gamma_array: 1st column, effective selection; 2nd column, proportion of sites
        
        gx = @(x,gamma,mu_site,L)2*mu_site*L*(1-exp(-1*gamma*(1-x)))/((1-exp(-1*gamma))*x*(1-x));

        function [total_SNPs,g] = m(gamma,mu_site,L,Npop)
            total_SNPs = 0;
            g = zeros((2*Npop-1),1);
            for j = 1:(2*Npop-1)
                freq = j/(2*Npop);
                if(gamma ~= 0) 
                    g(j) = gx(freq,gamma,mu_site,L); 
                else
                    g(j) = 2*mu_site*L/freq; 
                end
                total_SNPs = total_SNPs + g(j);
            end
        end

        function [prob] = prob(k,rho)
            if(length(P_samp) > length(rho))
                prob = dot(rho,P_samp((k+1),(2:length(P_samp))));
            else
                prob = dot(rho,P_samp((k+1),:));
            end

        end

        function [G,m_samp] = true_sample_G()
            m_mix = 0;
            g_mix = 0;

            for i = 1:length(gamma_array(:,1))
                [m_1,g_1] = m(gamma_array(i,1),mu_persite,(Num_sites*gamma_array(i,2)),N_pop);
                m_mix = m_mix + m_1;
                g_mix = g_mix + g_1;
            end

            m_samp = 0;
            rho_mix = g_mix/m_mix;
            G = zeros((N_samp-1),1);

            for k = 1:(N_samp-1)
                mix_prob = prob(k,rho_mix);
                %G(k) = round(mix_prob*m_mix);
                G(k) = mix_prob*m_mix;
                m_samp = m_samp + G(k);
            end       
        end

        [G,m_samp] = true_sample_G();   
    end

   function [loglambda,alpha,g_mix_sel_fold,g_mix_neu_fold] = total_likelihood(my_selected_sfs,my_neutral_sfs,P_samp,theta_site,p_array,gamma_array1,boundary_array,alpha,lethal_perc,zero_class,print_debug)
        p_array = abs(p_array);
        alpha = abs(alpha);
        theta_site = abs(theta_site);
        
        if(length(gamma_array1) == 1)
           p_array = [];
        end
        if((sum(p_array) > (1-lethal_perc)) || (lethal_perc > 1) || (lethal_perc < 0))
            loglambda = -1/0;
            return;
        end
        
        bool = 0;
        for index1 = 1:length(gamma_array1)
            if(boundary_array(index1,2) == boundary_array(index1,1))
                gamma_array1(index1) = boundary_array(index1,2);
            else
                bool = bool + (gamma_array1(index1) < boundary_array(index1,2)) + (gamma_array1(index1) > boundary_array(index1,1));
            end
        end
        if(bool>0)
           % bool,gamma_array1
            loglambda = -1/0;
            return;
        end 
        
        Num_selected_sites = sum(my_selected_sfs);
        Num_neutral_sites = sum(my_neutral_sfs);
        N_samp = 2*(length(my_selected_sfs)-1);
        mu_persite = theta_site/(4*N_pop);
        g_mix_neu = theor_SFS([0, 1.0],mu_persite,Num_neutral_sites,N_pop,N_samp,P_samp);

        %fold spectrum
        g_mix_neu_fold = zeros((length(my_neutral_sfs)-1),1);
        g_mix_neu_fold(1) = 1*(g_mix_neu(1) + g_mix_neu(N_samp - 1));
        for k = 2:(length(my_neutral_sfs)-1)
            g_mix_neu_fold(k) = alpha((k-1))*(g_mix_neu(k) + g_mix_neu(N_samp - k));
            if(k == N_samp - k) 
                g_mix_neu_fold(k) = alpha((k-1))*g_mix_neu(k); 
            end
        end    
        
        if(zero_class)
            %can't use m_samp because alphas, so re-sum number of neutral SNPs
            g0_mix_neu = Num_neutral_sites - sum(g_mix_neu_fold);
            g_mix_neu_fold = [g0_mix_neu; g_mix_neu_fold];
            rho_mix_neu_fold = g_mix_neu_fold/Num_neutral_sites;
        else
            rho_mix_neu = g_mix_neu_fold/sum(g_mix_neu_fold);
            rho_mix_neu_fold = [0;rho_mix_neu];
        end
             
        gamma_array_comb = [gamma_array1;[p_array,(1-lethal_perc)-sum(p_array)]]'; 
        g_mix_sel = theor_SFS(gamma_array_comb,mu_persite,Num_selected_sites,N_pop,N_samp,P_samp);
     
        if(print_debug) 
           % g_1
           %alpha
           %neutral_sfs
          % g_mix_neu(2:length(g_mix_neu))
        end
        
        g_mix_sel_fold = zeros((length(my_neutral_sfs)-1),1);
        g_mix_sel_fold(1) = 1*(g_mix_sel(1) + g_mix_sel(N_samp - 1));
        for k = 2:(length(my_neutral_sfs)-1)
            g_mix_sel_fold(k) = alpha((k-1))*(g_mix_sel(k) + g_mix_sel(N_samp - k));
            if(k == N_samp - k) 
                g_mix_sel_fold(k) = alpha((k-1))*g_mix_sel(k); 
            end
        end 
            
        if(zero_class)
            %can't use m_samp because alphas, so re-sum number of selected SNPs
            g0_mix_sel = Num_selected_sites - sum(g_mix_sel_fold);
            g_mix_sel_fold = [g0_mix_sel; g_mix_sel_fold];
            rho_mix_sel_fold = g_mix_sel_fold/Num_selected_sites;
        else
            rho_mix_sel = g_mix_sel_fold/sum(g_mix_sel_fold);
            rho_mix_sel_fold = [0;rho_mix_sel];
        end
        
        if(print_debug)
          %  g_mix_sel_fold
          %  g_mix_neu_fold
          %  m_mix_sel
          %  m_mix_neu
          %  wtf = Num_selected_sites - 3.5*power(10,5)
          %  wtf = Num_neutral_sites - 3.5*power(10,5)
          %  N_samp
          %  sum(selected_sfs(2:length(selected_sfs)))
          %  sum(neutral_sfs(2:length(neutral_sfs)))
        end
        
        loglambda_neu = 0;
        loglambda_sel = 0;
        
        if(zero_class)
            loglambda_neu = my_neutral_sfs(1)*log(rho_mix_neu_fold(1));
            loglambda_sel = my_selected_sfs(1)*log(rho_mix_sel_fold(1)); 
        end
           
        
        for k = 2:length(my_neutral_sfs)
            loglambda_neu = loglambda_neu + my_neutral_sfs(k)*log(rho_mix_neu_fold(k));
        end    
            
        for k = 2:length(my_selected_sfs)
            loglambda_sel = loglambda_sel + my_selected_sfs(k)*log(rho_mix_sel_fold(k));
        end
        
        loglambda = loglambda_neu + loglambda_sel;
   end

    function [loglambda] = re_param_neutral_alpha(x)
        theta_site = x(7);
        p_array = [];
        
        alpha(1) = x(1);
        alpha(2:3) = x(2);
        alpha(4:7) = x(3);
        alpha(8:15) = x(4);
        alpha(16:  floor((length(neutral_sfs)-18)/2) ) = x(5);
        alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = x(6);    

        gamma_array_1 = 0;
        b_array = [0,0];       
        loglambda = -1*total_likelihood(neutral_sfs,neutral_sfs,P__samp,theta_site,p_array,gamma_array_1,b_array,alpha,0,zero_class,false);
    end
    
    function [loglambda] = re_param(x, fixed_alpha)
        theta_site = x(7);
        p_array = x(8: (8+length(initial_p_array)-1));
        
        if(free_alpha && ~neutral_alpha)
            alpha(1) = x(1);
            alpha(2:3) = x(2);
            alpha(4:7) = x(3);
            alpha(8:15) = x(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = x(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = x(6);    
        else
            alpha(1) = fixed_alpha(1);
            alpha(2:3) = fixed_alpha(2);
            alpha(4:7) = fixed_alpha(3);
            alpha(8:15) = fixed_alpha(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = fixed_alpha(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = fixed_alpha(6);
        end

        gamma_array_1 = x( (8+length(p_array)):length(x) );

        loglambda = -1*total_likelihood(selected_sfs,neutral_sfs,P__samp,theta_site,p_array,gamma_array_1,boundary_array,alpha,0,zero_class,false);
    end

    function [loglambda] = re_param_neu(x, fixed_alpha)
        theta_site = x(7);
        p_array = [];
        
        if(free_alpha && ~neutral_alpha)
            alpha(1) = x(1);
            alpha(2:3) = x(2);
            alpha(4:7) = x(3);
            alpha(8:15) = x(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = x(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = x(6);    
        else
            alpha(1) = fixed_alpha(1);
            alpha(2:3) = fixed_alpha(2);
            alpha(4:7) = fixed_alpha(3);
            alpha(8:15) = fixed_alpha(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = fixed_alpha(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = fixed_alpha(6);
        end

        gamma_array_1 = 0;
        b_array = [0,0];
        loglambda = -1*total_likelihood(selected_sfs,neutral_sfs,P__samp,theta_site,p_array,gamma_array_1,b_array,alpha,0,zero_class,false);
    end

    function [loglambda] = re_param_neu_lethal(x, fixed_alpha)
        theta_site = x(7);
        p_array = 1.0-x(length(x));
        
        if(free_alpha && ~neutral_alpha)
            alpha(1) = x(1);
            alpha(2:3) = x(2);
            alpha(4:7) = x(3);
            alpha(8:15) = x(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = x(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = x(6);    
        else
            alpha(1) = fixed_alpha(1);
            alpha(2:3) = fixed_alpha(2);
            alpha(4:7) = fixed_alpha(3);
            alpha(8:15) = fixed_alpha(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = fixed_alpha(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = fixed_alpha(6);
        end
        
        gamma_array_1 = 0;
        b_array = [0,0];
        lethal_perc = x(length(x));
        loglambda = -1*total_likelihood(selected_sfs,neutral_sfs,P__samp,theta_site,p_array,gamma_array_1,b_array,alpha,lethal_perc,zero_class,false);
    end

    function [loglambda] = re_param_lethal(x, fixed_alpha)
        theta_site = x(7);
        p_array = x(8: (8+length(initial_p_array)-1));
        
        if(free_alpha && ~neutral_alpha)
            alpha(1) = x(1);
            alpha(2:3) = x(2);
            alpha(4:7) = x(3);
            alpha(8:15) = x(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = x(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = x(6);    
        else
            alpha(1) = fixed_alpha(1);
            alpha(2:3) = fixed_alpha(2);
            alpha(4:7) = fixed_alpha(3);
            alpha(8:15) = fixed_alpha(4);
            alpha(16:  floor((length(neutral_sfs)-18)/2) ) = fixed_alpha(5);
            alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = fixed_alpha(6);
        end

        gamma_array_1 = x( (8+length(p_array)):(length(x)-1) );
        lethal_perc = x(length(x));
        
        loglambda = -1*total_likelihood(selected_sfs,neutral_sfs,P__samp,theta_site,p_array,gamma_array_1,boundary_array,alpha,lethal_perc,zero_class,false);
    end

    x_neu_alpha = [ones(1,6),initial_theta_site];
    old_set = optimset('fminsearch');
    new_set = optimset(old_set, 'MaxFunEvals', 100000*(length(x_neu_alpha)+length(initial_p_array)+length(initial_gamma_array)),'MaxIter', 100000*(length(x_neu_alpha)+length(initial_p_array)+length(initial_gamma_array)));
    init_alpha = ones(1,6);
    
    if(neutral_alpha && free_alpha)
        x_max_neu_alpha = fminsearch(@re_param_neutral_alpha,x_neu_alpha,new_set);
        init_alpha = x_max_neu_alpha(1:6);
    end

    x = [ones(1,6),initial_theta_site,initial_p_array,initial_gamma_array]; 

    x_max = fminsearch(@(x) re_param(x,init_alpha),x,new_set);
 
    if(free_alpha && ~neutral_alpha)
        max_alpha(1) = x_max(1); 
        max_alpha(2:3) = x_max(2);
        max_alpha(4:7) = x_max(3);
        max_alpha(8:15) = x_max(4);
        max_alpha(16:  floor((length(neutral_sfs)-18)/2) ) = x_max(5);
        max_alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = x_max(6);  
    else
        max_alpha(1) = init_alpha(1); 
        max_alpha(2:3) = init_alpha(2);
        max_alpha(4:7) = init_alpha(3);
        max_alpha(8:15) = init_alpha(4);
        max_alpha(16:  floor((length(neutral_sfs)-18)/2) ) = init_alpha(5);
        max_alpha( floor((length(neutral_sfs)-17) / 2) : (length(neutral_sfs)-2) ) = init_alpha(6); 
    end
    
    max_theta_site = abs(x_max(7));   
    max_p_array = abs(x_max(8: (8+length(p_array)-1) ));   
    max_gamma = x_max( (8+length(p_array)):length(x) );  
    [max_likelihood,max_alpha,g_mix_sel_fold,g_mix_neu_fold] = total_likelihood(selected_sfs,neutral_sfs,P__samp,max_theta_site,max_p_array,max_gamma,boundary_array,max_alpha,0,zero_class,false);
    
    x_neu = [ones(1,6),initial_theta_site];
    [x_max_neu,all_neutral_lik] = fminsearch(@(x_neu) re_param_neu(x_neu,init_alpha),x_neu,new_set);    
    all_neutral_lik = -1*all_neutral_lik;
    max_theta_site_neu = x_max_neu(7);
    
    if(zero_class)
        x_neu_lethal = [ones(1,6),initial_theta_site,initial_lethal_array(1)]; 
        [x_max_neu_lethal,max_neutral_lethal_lik] = fminsearch(@(x_neu_lethal) re_param_neu_lethal(x_neu_lethal,init_alpha),x_neu_lethal,new_set);    
        max_neutral_lethal_lik = -1*max_neutral_lethal_lik;
        max_neu_lethal_perc = x_max_neu_lethal(length(x_max_neu_lethal));
        max_theta_site_neu_lethal = x_max_neu_lethal(7);
        
        x_lethal = [ones(1,6),initial_theta_site,initial_p_array,initial_gamma_array,initial_lethal_array(2)]; 
        [x_max_lethal,max_lethal_lik] = fminsearch(@(x_lethal) re_param_lethal(x_lethal,init_alpha),x_lethal,new_set);
        max_lethal_perc = x_max_lethal(length(x_max_lethal));
        max_lethal_lik = -1*max_lethal_lik;
        max_theta_site_lethal = x_max_lethal(7);
        max_lethal_p_array = abs(x_max_lethal(8: (8+length(initial_p_array)-1)));
        max_lethal_gamma = x_max_lethal( (8+length(initial_p_array)):(length(x_max_lethal)-1) );
    end

    fileID = fopen(file_name, 'A');
    fprintf(fileID, '%d_selection_classes\t', length(initial_gamma_array));
   
    if(zero_class)
        fprintf(fileID, '%s\t', 'zero_class');
    else
        fprintf(fileID, '%s\t', 'shape_only');
    end

    if(free_alpha && ~neutral_alpha)
        fprintf(fileID, '%s\t', 'freeAlpha6Bins');
    elseif(neutral_alpha && free_alpha)
        fprintf(fileID, '%s\t', 'neutAlpha6Bins');
    else
        fprintf(fileID, '%s\t', 'fixedAlpha_one');
    end
    
    fprintf(fileID, '(');
    for index = 1:(length(initial_gamma_array)-1)
        fprintf(fileID, '%s,', num2str(initial_gamma_array(index)));
    end
    fprintf(fileID, '%s;', num2str(initial_gamma_array(length(initial_gamma_array))));
    
    for index = 1:(length(max_gamma)-1)
        fprintf(fileID, '%s,', mat2str(boundary_array(index,:)));
    end
    fprintf(fileID, '%s;', mat2str(boundary_array(length(max_gamma),:)));
    
    for index = 1:(length(initial_p_array))
        fprintf(fileID, '%s,', num2str(initial_p_array(index)));
    end
    fprintf(fileID, '%s;',num2str((1-sum(initial_p_array))));
    
    fprintf(fileID, '%6.4f;', initial_theta_site);
    fprintf(fileID, '%s,%s,%s', num2str(initial_lethal_array(1)), num2str(initial_lethal_array(2)), num2str(N_pop));
    fprintf(fileID, ')\t');
    
    fprintf(fileID, '(');
    for index = 1:(length(max_gamma)-1)
        if(boundary_array(index,2) == boundary_array(index,1)) 
            max_gamma(index) = boundary_array(index,2);
        end
        fprintf(fileID, '%s %6.4f, ', num2str(max_gamma(index)), max_p_array(index));
    end  
    if(boundary_array(length(max_gamma),2) == boundary_array(length(max_gamma),1)) 
        max_gamma(length(max_gamma)) = boundary_array(length(max_gamma),2);
    end
    pmax = 1-sum(max_p_array);
    if(length(max_gamma) == 1) 
        pmax = 1;
    end
    if(zero_class)
        fprintf(fileID, '%s %6.4f, theta %6.4f): ', num2str(max_gamma(length(max_gamma))), pmax, max_theta_site);
    else
        fprintf(fileID, '%s %6.4f): ', num2str(max_gamma(length(max_gamma))), pmax);
    end
    fprintf(fileID,'%6.4f\t', max_likelihood);
    
    if(zero_class)
        fprintf(fileID,'(neutral 1.0, theta %6.4f): %6.4f\t', max_theta_site_neu, all_neutral_lik);
    else
        fprintf(fileID,'(neutral 1.0): %6.4f\t', all_neutral_lik);
    end
    
    if(zero_class)
        fprintf(fileID,'(neutral %6.4f, lethal %6.4f, theta %6.4f): %6.4f\t', (1-max_neu_lethal_perc), max_neu_lethal_perc, max_theta_site_neu_lethal, max_neutral_lethal_lik);
        
        fprintf(fileID, '(');    
        for index = 1:(length(max_lethal_gamma)-1)
            if(boundary_array(index,2) == boundary_array(index,1)) 
                max_lethal_gamma(index) = boundary_array(index,2);
            end
            fprintf(fileID, '%s %6.4f, ', num2str(max_lethal_gamma(index)), max_lethal_p_array(index));
        end
        if(boundary_array(length(max_lethal_gamma),2) == boundary_array(length(max_lethal_gamma),1)) 
            max_lethal_gamma(length(max_lethal_gamma)) = boundary_array(length(max_lethal_gamma),2);
        end
        pmax = 1-sum(max_lethal_p_array)-max_lethal_perc;
        if(length(max_gamma) == 1) 
            pmax = 1-max_lethal_perc;
        end
        fprintf(fileID,'%s %6.4f, lethal %6.4f, theta %6.4f): %6.4f', num2str(max_lethal_gamma(length(max_lethal_gamma))), pmax, max_lethal_perc, max_theta_site_lethal, max_lethal_lik);
    end
    fprintf(fileID, '\n');
    fclose(fileID);

end