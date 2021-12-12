function longterm_projecct4

        weights_from_S_to_L4 = [];
        weights_from_D_to_L4 = [];
        weights_from_SP_to_L4 = [];

        weight_S_to_SP = 0.2;
        weight_D_to_SP = 0.2;

        weight_SP_to_L4 = 0.11;
        weight_S_to_L4 = 0.02;
        weight_D_to_L4 = 0.02;

        a_LTP = 0.015;
        tau_LTP = 13;

        a_LTD = 0.021;
        tau_LTD = 20;


        for p=1:180 % for 9 minutes for now
            stimulus_decider = randi([1,100]);
            spike_train_of_thalamus_S = [];
            spike_train_of_thalamus_D = [];

            voltage_for_sp_for_single_stimulus = [];

            if stimulus_decider >= 90
                % Deviant stimulus
                spike_for_S = binornd(1, 2.5/1000, 1,1000);
                spike_for_D = binornd(1, 10/1000, 1,1000);
            else 
                % Standard stimulus
                spike_for_S = binornd(1, 10/1000, 1,1000);
                spike_for_D = binornd(1, 2.5/1000, 1,1000);
            end

            % 50 ms stimulus
            k_t1 = get_kernel2(spike_for_S);
            g_t1 = conv(k_t1, spike_for_S);
            g_t1 = g_t1(1,1:50);
            g_t1 = shift_by_2(g_t1);

            k_t2 = get_kernel2(spike_for_D);
            g_t2 = conv(k_t2, spike_for_D);
            g_t2 = g_t2(1,1:50);
            g_t2 = shift_by_2(g_t2);


            spike_for_S = spike_for_S(1,1:50);
            spike_for_S = shift_by_2(spike_for_S);
            
            spike_for_D = spike_for_D(1, 1:50);
            spike_for_D = shift_by_2(spike_for_D);
        
        
        spike_train_of_thalamus_S = [spike_train_of_thalamus_S, spike_for_S];
        spike_train_of_thalamus_D = [spike_train_of_thalamus_D, spike_for_D];

        voltage_when_stimulus = g_t1.*(weight_S_to_SP*spike_for_S) + g_t2.*(weight_D_to_SP*spike_for_D);
        voltage_for_sp_for_single_stimulus = [voltage_for_sp_for_single_stimulus, voltage_when_stimulus];
   

        % 50 ms no stimulus
        spike_for_S = binornd(1, 0.5/1000, 1,1000);
        spike_for_D = binornd(1, 0.5/1000, 1,1000);
        
        

        k_t1 = get_kernel(spike_for_S);
        g_t1 = conv(k_t1, spike_for_S);
        g_t1 = g_t1(1,1:250);
        g_t1 = shift_by_2(g_t1);

        k_t2 = get_kernel(spike_for_D);
        g_t2 = conv(k_t2, spike_for_D);
        g_t2 = g_t2(1,1:250);
        g_t2 = shift_by_2(g_t2);

        spike_for_S = spike_for_S(1,1:250);
        spike_for_S = shift_by_2(spike_for_S);
        
        spike_for_D = spike_for_D(1, 1:250);
        spike_for_D = shift_by_2(spike_for_D);

        spike_train_of_thalamus_S = [spike_train_of_thalamus_S, spike_for_S];
        spike_train_of_thalamus_D = [spike_train_of_thalamus_D, spike_for_D];


        voltage_when_no_stimulus = g_t1.*(weight_S_to_SP*spike_for_S) + g_t2.*(weight_D_to_SP*spike_for_D);
        voltage_for_sp_for_single_stimulus= [ voltage_for_sp_for_single_stimulus, voltage_when_no_stimulus];
        
         % decrease 10% for every voltage to account for leak 
        for i=1:300
            voltage_for_sp_for_single_stimulus(1,i) = voltage_for_sp_for_single_stimulus(1,i) - (0.1*voltage_for_sp_for_single_stimulus(1,i));
        end

    [voltage_for_sp_total, spike_train_for_SP] = decrease_voltage_for_20ms_after_spike(voltage_for_sp_for_single_stimulus);



        % L4 voltage
        % from S, from D, from SP
        k_t3 = get_kernel2(spike_train_of_thalamus_S);
        g_t3 = conv(k_t3, spike_train_of_thalamus_S);
        g_t3 = g_t3(1,1:300);
        g_t3 = shift_by_2(g_t3);

        k_t4 = get_kernel(spike_train_of_thalamus_D);
        g_t4 = conv(k_t4, spike_train_of_thalamus_D);
        g_t4 = g_t4(1,1:300);
        g_t4 = shift_by_2(g_t4);

        
        k_t5 = get_kernel(spike_train_for_SP);
        g_t5 = conv(k_t5, spike_train_for_SP);
        g_t5 = g_t5(1,1:300);
        g_t5 = shift_by_2(g_t5);

        spike_train_of_thalamus_S = spike_train_of_thalamus_S(1,1:300);
        spike_train_of_thalamus_S = shift_by_2(spike_train_of_thalamus_S);

        spike_train_of_thalamus_D = spike_train_of_thalamus_D(1,1:300);
        spike_train_of_thalamus_D = shift_by_2(spike_train_of_thalamus_D);

    
        voltage_for_L4 = g_t3.*(weight_S_to_L4*spike_train_of_thalamus_S) + g_t4.*(weight_D_to_L4*spike_train_of_thalamus_D) + g_t5.*(weight_SP_to_L4*spike_train_for_SP);
        
        [voltage_for_L4, spike_train_for_L4] = decrease_voltage_for_20ms_after_spike(voltage_for_L4);

        % TODO
        % 
        % - change weights at the end for the next stimuli

        weights_from_S_to_L4 = [weights_from_S_to_L4, weight_S_to_L4];
        weights_from_D_to_L4 = [weights_from_D_to_L4, weight_D_to_L4];
        weights_from_SP_to_L4 = [weights_from_SP_to_L4, weight_SP_to_L4];

        % spike_train_for_L4
        % spike_train_for_SP
        % spike_train_of_thalamus_S
        % spike_train_of_thalamus_D

        
        % what hapens to connections from SP to L4\
        % to know weights effect
        voltage_L4_due_to_SP = g_t5.*(weight_SP_to_L4*spike_train_for_SP);
        [voltage_L4_due_to_SP, spike_train_L4_due_to_SP] =  decrease_voltage_for_20ms_after_spike(voltage_L4_due_to_SP);
        effective_presyn_spikes_in_SP = 0;
        for i=1:300
            if spike_train_for_SP(1,i) == 1
                for j=i+1:i+10 % checking if any post synaptic spike caused within 10ms
                    if j < 300
                        if spike_train_L4_due_to_SP(1,j) == 1
                            effective_presyn_spikes_in_SP = effective_presyn_spikes_in_SP + 1;
                            break
                        end
                    end
                end
            end
        end

        fprintf("\n eff is %d \n", effective_presyn_spikes_in_SP);
        
        if effective_presyn_spikes_in_SP >= 150
            weight_SP_to_L4 = weight_SP_to_L4 * ( 1 + (a_LTP * exp(-(10)/tau_LTP)) );
        else
            weight_SP_to_L4 = weight_SP_to_L4 * ( 1 + (-a_LTD * exp(-(10)/tau_LTD)) );
        end

        if weight_SP_to_L4 <= 0.001
            weight_SP_to_L4 = 0.0001;
        end

        if weight_SP_to_L4 >= 0.11
            weight_SP_to_L4 = 0.11;
        end




        % read figure
        % what happens to connections from S, D to L4


        end % 1800 Stimulus - 9 m
end

function  [new_voltage_values spike_train] = decrease_voltage_for_20ms_after_spike(voltage_values)
    voltage_after_decreasing = zeros(1, 300);
    nearest_spike_time = 0;

    actual_spike = zeros(1, 300);

    % variables for v_delta
    beta = 5; tau = 2;

    for i=1:300
        if voltage_values(1,i) < 0.05
            voltage_after_decreasing(1,i) = voltage_values(1,i);
            continue
        end

        if voltage_values(1,i) >= 0.05
            nearest_spike_time = i;
            for j=i:i+19
                if j < 300
                    voltage_after_decreasing(1,j) = voltage_values(1,j)*(1 - ( beta * exp(-(j-nearest_spike_time)/tau) ))  ; 
                    actual_spike(1,j) = 1;
                end
            end
        end

        i = j;
        
    end


    new_voltage_values = voltage_after_decreasing;
    spike_train = actual_spike;
end

function shifted_arr = shift_by_2(arr)
    len = length(arr);
    new_arr = zeros(1, len);
    for i=3:len
        new_arr(1,i) = arr(1, i-2); 
    end

    shifted_arr = new_arr;
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