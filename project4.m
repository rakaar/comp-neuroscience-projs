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

    % spikes from S and D for calculating Spike of L4
    spike_train_of_thalamus_S = [];
    spike_train_of_thalamus_D = [];
    for i=1:15
        spike_for_S = binornd(1, 10/1000, 1,1000);
        spike_for_D = binornd(1, 2.5/1000, 1,1000);
        if i == 8
            spike_for_S = binornd(1, 2.5/1000, 1,1000);
            spike_for_D = binornd(1, 10/1000, 1,1000);
        end
       
        spike_train_of_thalamus_S = [spike_train_of_thalamus_S, spike_for_S];
        spike_train_of_thalamus_D = [spike_train_of_thalamus_D, spike_for_D];

        voltage_for_sp_for_single_stimulus = zeros(1,300);
        voltage_for_sp_for_single_stimulus(1,1) = -0.070;
        
        
        k_t1 = get_kernel(spike_for_S);
        g_t1 = conv(k_t1, spike_for_S);
        g_t1 = g_t1(1,1:250);

        k_t2 = get_kernel(spike_for_D);
        g_t2 = conv(k_t2, spike_for_D);
        g_t2 = g_t2(1,1:250);

        spike_for_S = spike_for_S(1,1:250);
        spike_for_D = spike_for_D(1, 1:250);
        
        spike_train_of_thalamus_S = [spike_train_of_thalamus_S, spike_for_S];
        spike_train_of_thalamus_D = [spike_train_of_thalamus_D, spike_for_D];

        voltage_when_stimulus = g_t1.*(weight_S_to_SP*spike_for_S) + g_t2.*(weight_D_to_SP*spike_for_D);
        
        for ind=2:251
            voltage_for_sp_for_single_stimulus(1, ind) = voltage_when_stimulus(1, ind-1);
        end

        % 50 ms no stimulus
        spike_for_S = binornd(1, 0.5/1000, 1,1000);
        spike_for_D = binornd(1, 0.5/1000, 1,1000);
        
        spike_for_S = spike_for_S(1,1:50);
        spike_for_D = spike_for_S(1,1:50);

        k_t1 = get_kernel(spike_for_S);
        g_t1 = conv(k_t1, spike_for_S);
        g_t1 = g_t1(1,1:50);

        k_t2 = get_kernel(spike_for_D);
        g_t2 = conv(k_t2, spike_for_D);
        g_t2 = g_t2(1,1:50);

        voltage_when_no_stimulus = g_t1.*(weight_S_to_SP*spike_for_S) + g_t2.*(weight_D_to_SP*spike_for_D);
        for ind=252:300
            voltage_for_sp_for_single_stimulus(1, ind) = voltage_when_no_stimulus(1, ind-251);
        end

        voltage_for_sp_total = [voltage_for_sp_total, voltage_for_sp_for_single_stimulus];

    end 

    figure(3)
        stem(voltage_for_sp_total);
    grid
   
    % decrease 10% for every voltage to account for leak 
    for i=1:1800
        voltage_for_sp_total(1,i) = voltage_for_sp_total(1,i) - (0.1*voltage_for_sp_total(1,i));
    end
    figure(5)
        stem(voltage_for_sp_total);
    grid

    % from  voltage to spike train for SP 
    voltage_for_sp_total = decrease_voltage_for_20ms_after_spike(voltage_for_sp_total);

    spike_train_for_SP = zeros(1,1800); 
    for i=1:1800
        if voltage_for_sp_total(1,i) >= 0.05
            spike_train_for_SP(1,i) = 1;
        end
    end

    % ------- spike train for SP neuron ------
    figure(8)
        stem(spike_train_for_SP);
    grid

    % L4 voltage
    % from S, from D, from SP
    % voltage_for_L4 = ;
    k_t3 = get_kernel2(spike_train_of_thalamus_S);
    g_t3 = conv(k_t3, spike_train_of_thalamus_S);
    g_t3 = g_t3(1,1:1800);

    k_t4 = get_kernel(spike_train_of_thalamus_D);
    g_t4 = conv(k_t4, spike_train_of_thalamus_D);
    g_t4 = g_t4(1,1:1800);
    
    k_t5 = get_kernel(spike_train_for_SP);
    g_t5 = conv(k_t5, spike_train_for_SP);
    g_t5 = g_t5(1,1:1800);

    spike_train_of_thalamus_S = spike_train_of_thalamus_S(1,1:1800);
    spike_train_of_thalamus_D = spike_train_of_thalamus_D(1,1:1800);
    % fprintf("\n length(g_t3) %d \n", length(g_t3));
    % fprintf("\n length(g_t4) %d \n", length(g_t4));
    % fprintf("\n length(g_t5) %d \n", length(g_t5));

    % fprintf("\n length(spike_train_of_thalamus_S) %d \n", length(spike_train_of_thalamus_S));
    % fprintf("\n length(spike_train_of_thalamus_D) %d \n", length(spike_train_of_thalamus_D));
    % fprintf("\n length(spike_train_for_SP) %d \n", length(spike_train_for_SP));


    voltage_for_L4 = g_t3.*(weight_S_to_L4*spike_train_of_thalamus_S) + g_t4.*(weight_D_to_L4*spike_train_of_thalamus_D) + g_t5.*(weight_SP_to_L4*spike_train_for_SP);
    
    % decrease 10% for every voltage to account for leak 
    for i=1:1800
        voltage_for_L4(1,i) = voltage_for_L4(1,i) - (0.1*voltage_for_L4(1,i));
    end
    figure(52)
        stem(voltage_for_L4);
    grid

    % from  voltage to spike train for SP 
    voltage_for_L4 = decrease_voltage_for_20ms_after_spike(voltage_for_L4);

    figure(54)
        stem(voltage_for_L4);
    grid
    spike_train_for_L4 = zeros(1,1800); 
    for i=1:1800
        if voltage_for_L4(1,i) >= 0.05
            spike_train_for_L4(1,i) = 1;
        end
    end

    figure(55)
        stem(spike_train_for_L4);
    grid

end

function  new_voltage_values = decrease_voltage_for_20ms_after_spike(voltage_values)
    spike_after_decreasing_voltage = zeros(1, 1800);
    nearest_spike_time = 0;

    % variables for v_delta
    beta = 5; tau = 2;

    for i=1:1800
        if voltage_values(1,i) < 0.05
            spike_after_decreasing_voltage(1,i) = voltage_values(1,i);
            continue
        end

        if voltage_values(1,i) >= 0.05
            nearest_spike_time = i;
            for j=i+1:i+20
                if j < 1800
                    spike_after_decreasing_voltage(1,j) = voltage_values(1,j)*(1 - ( beta * exp(-(j-nearest_spike_time)/tau) ))  ; 
                end
            end
        end

        i = j;
        
    end


    new_voltage_values = spike_after_decreasing_voltage;
end


function kernel2 = get_kernel2(spike_train)
    neuron_kernel = zeros(1,length(spike_train));
    nearest_spike_time = 0;
    for i=1:length(spike_train)
        if spike_train(1,i) == 1
            nearest_spike_time = i;
        end

        neuron_kernel(1,i) = exp(- (i - nearest_spike_time)/10);
    end

    kernel2 = neuron_kernel;
end


function kernel = get_kernel(spike_train)
    neuron_kernel = zeros(1,1000);
    nearest_spike_time = 0;
    for i=1:length(spike_train)
        if spike_train(1,i) == 1
            nearest_spike_time = i;
        end
        neuron_kernel(1,i) = exp(- (i - nearest_spike_time)/10);
    end

    kernel = neuron_kernel;
end

