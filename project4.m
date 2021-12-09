function project

    %% non homo - testing of non homogenorus
    % a poisson distribution with just 0s and 1s that is non homogeneous = binomial distrubutiion with varying probability
    random_small_prob = randi([40,50], 100,1000)/1000;
    nonhomopp = binornd(1, random_small_prob, 100,1000);

    disp(max(nonhomopp, [], 'all'))
    
    % generate a homogenous poison process with different repitions- reps = [10, 20,40,80]
    % TODO - RMSE of exactly what?? PSTH is 1x1000, the reps  are 100x1000, how rmse?
    
    PSTH = zeros(1,1000);
    for j=1:1000
        for i=1:100
            PSTH(1,j) = PSTH(1,j) + nonhomopp(i,j); 
        end
        PSTH(1,j) = PSTH(1,j)/(100*0.001);
    end

    figure(1)
        plot(linspace(1,1000,1000), PSTH);
    grid

    % -------- Short Term Plasticity ------

    %% generate spike times for S and D neurons
    no_stimulus_for_S = binornd(1, 0.5/1000, 1,1000); 
    S_stimulus_for_thalamus_S = binornd(1, 10/1000, 1,1000); 
    D_stimulus_for_thalamus_S = binornd(1, 2.5/1000, 1,1000); 
    
    
    no_stimulus_for_D = binornd(1, 0.5/1000, 1,1000);  
    S_stimulus_for_thalamus_D = binornd(1, 2.5/1000, 1,1000); 
    D_stimulus_for_thalamus_D = binornd(1, 10/1000, 1,1000); 

    weight_S_to_SP = 0.2;
    weight_D_to_SP = 0.2;
    weight_SP_to_L4 = 0.11;
    weight_S_to_L4 = 0.02;
    weight_D_to_L4 = 0.02;


    %% generate voltage for SP neuron
    voltage_for_sp_total = []; %- length 1800 -  15 stimuli * (250 ms stimulus played +  50 ms gap)
    for i=1:15
        voltage_for_sp_for_single_stimulus = zeros(1,300);
        voltage_for_sp_for_single_stimulus(1,1) = -0.070;
       
        spike_for_S = binornd(1, 10/1000, 1,1000);
        spike_for_D = binornd(1, 2.5/1000, 1,1000);
        if i == 8
            spike_for_S = binornd(1, 2.5/1000, 1,1000);
            spike_for_D = binornd(1, 10/1000, 1,1000);
        end
       
        k_t1 = get_kernel(spike_for_S);
        g_t1 = conv(k_t1, spike_for_S);
        g_t1 = g_t1(1,1:300);

        k_t2 = get_kernel(spike_for_D);
        g_t2 = conv(k_t2, spike_for_D);
        g_t2 = g_t2(1,1:300);

        spike_for_S = spike_for_S(1,1:300);
        spike_for_D = spike_for_D(1, 1:300);

        voltage_from_t_is_2 = g_t1.*(weight_S_to_SP*spike_for_S) + g_t2.*(weight_D_to_SP*spike_for_D);
        
        for ind=2:300
            voltage_for_sp_for_single_stimulus(1, ind) = voltage_from_t_is_2(1, ind-1);
        end

        
        voltage_for_sp_total = [voltage_for_sp_total, voltage_for_sp_for_single_stimulus];

    end 

    figure(3)
        stem(voltage_for_sp_total);
    grid
   

    % from  voltage to spike train for SP 


    % L4 voltage

    % L4 spike train

end

function kernel = get_kernel(spike_train)
    neuron_kernel = zeros(1,1000);
    nearest_spike_time = 0;
    for i=1:1000
        if spike_train(1,i) == 1
            nearest_spike_time = i;
        end
        neuron_kernel(1,i) = exp(- (i - nearest_spike_time)/10);
    end

    kernel = neuron_kernel;
end

