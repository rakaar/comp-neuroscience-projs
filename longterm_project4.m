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


        for p=1:50 % for 9 minutes for now
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
            k_t1 = shift_by_x(k_t1, 1);
            g_t1 = conv(k_t1, spike_for_S);
            g_t1 = g_t1(1,1:50);
            g_t1 = shift_by_x(g_t1,1);

            k_t2 = get_kernel2(spike_for_D);
            k_t2 = shift_by_x(k_t2, 1);
            g_t2 = conv(k_t2, spike_for_D);
            g_t2 = g_t2(1,1:50);
            g_t2 = shift_by_x(g_t2,1);


            spike_for_S = spike_for_S(1,1:50);
            spike_for_S = shift_by_x(spike_for_S,1);
            
            spike_for_D = spike_for_D(1, 1:50);
            spike_for_D = shift_by_x(spike_for_D,1);
        
        
        spike_train_of_thalamus_S = [spike_train_of_thalamus_S, spike_for_S];
        spike_train_of_thalamus_D = [spike_train_of_thalamus_D, spike_for_D];

        voltage_when_stimulus = g_t1.*(weight_S_to_SP*spike_for_S) + g_t2.*(weight_D_to_SP*spike_for_D);
        voltage_for_sp_for_single_stimulus = [voltage_for_sp_for_single_stimulus, voltage_when_stimulus];
   

        % 50 ms no stimulus
        spike_for_S = binornd(1, 0.5/1000, 1,1000);
        spike_for_D = binornd(1, 0.5/1000, 1,1000);
        
        

        k_t1 = get_kernel(spike_for_S);
        k_t1 = shift_by_x(k_t1, 1);
        g_t1 = conv(k_t1, spike_for_S);
        g_t1 = g_t1(1,1:250);
        g_t1 = shift_by_x(g_t1,1);

        k_t2 = get_kernel(spike_for_D);
        k_t2 = shift_by_x(k_t2, 1);
        g_t2 = conv(k_t2, spike_for_D);
        g_t2 = g_t2(1,1:250);
        g_t2 = shift_by_x(g_t2,1);

        spike_for_S = spike_for_S(1,1:250);
        spike_for_S = shift_by_x(spike_for_S,1);
        
        spike_for_D = spike_for_D(1, 1:250);
        spike_for_D = shift_by_x(spike_for_D,1);

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
        k_t3 = shift_by_x(k_t3, 1);
        g_t3 = conv(k_t3, spike_train_of_thalamus_S);
        g_t3 = g_t3(1,1:300);
        g_t3 = shift_by_x(g_t3,1);

        k_t4 = get_kernel(spike_train_of_thalamus_D);
        k_t4 = shift_by_x(k_t4, 1);
        g_t4 = conv(k_t4, spike_train_of_thalamus_D);
        g_t4 = g_t4(1,1:300);
        g_t4 = shift_by_x(g_t4,1);

        
        k_t5 = get_kernel(spike_train_for_SP);
        k_t5 = shift_by_x(k_t5, 1);
        g_t5 = conv(k_t5, spike_train_for_SP);
        g_t5 = g_t5(1,1:300);
        g_t5 = shift_by_x(g_t5,1);

        spike_train_of_thalamus_S = spike_train_of_thalamus_S(1,1:300);
        spike_train_of_thalamus_S = shift_by_x(spike_train_of_thalamus_S,1);

        spike_train_of_thalamus_D = spike_train_of_thalamus_D(1,1:300);
        spike_train_of_thalamus_D = shift_by_x(spike_train_of_thalamus_D,1);

    
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

        
        %-------- connections from SP to L4
        voltage_L4_due_to_SP = g_t5.*(weight_SP_to_L4*spike_train_for_SP);
        [voltage_L4_due_to_SP, spike_train_L4_due_to_SP] =  decrease_voltage_for_20ms_after_spike(voltage_L4_due_to_SP);
        arr_of_delta_ts = get_post_minus_pre_array(spike_train_for_SP, spike_train_L4_due_to_SP); 
        for i=1:length(arr_of_delta_ts)
            if arr_of_delta_ts(1,i) > 0
                weight_SP_to_L4 = weight_SP_to_L4 * ( 1 + (a_LTP * exp(-arr_of_delta_ts(1,i)/tau_LTP)) );
            else
                weight_SP_to_L4 = weight_SP_to_L4 * ( 1 + (-a_LTD * exp(arr_of_delta_ts(1,i)/tau_LTD)) );
            end
        end
        % make sure weight doesn't cross limits
        if weight_SP_to_L4 <= 0.001
            weight_SP_to_L4 = 0.0001;
        end
        if weight_SP_to_L4 >= 0.11
            weight_SP_to_L4 = 0.11;
        end


        %-------- connections from S to L4
        voltage_L4_due_to_S = g_t3.*(weight_S_to_L4*spike_train_of_thalamus_S);
        [voltage_L4_due_to_S, spike_train_L4_due_to_S] =  decrease_voltage_for_20ms_after_spike(voltage_L4_due_to_S);
        
        arr_of_delta_ts = get_post_minus_pre_array(spike_train_of_thalamus_S, spike_train_L4_due_to_S);

        disp("Arr coming?")
        disp(max(spike_train_of_thalamus_S));
        disp(max(spike_train_L4_due_to_S));
        disp(arr_of_delta_ts);
        for i=1:length(arr_of_delta_ts)
            if arr_of_delta_ts(1,i) > 0
                weight_S_to_L4 = weight_S_to_L4 * ( 1 + (a_LTP * exp(-arr_of_delta_ts(1,i)/tau_LTP)) );
            else
                weight_S_to_L4 = weight_S_to_L4 * ( 1 + (-a_LTD * exp(arr_of_delta_ts(1,i)/tau_LTD)) );
            end
        end

        
        fprintf("\n  weight s - l4 %f  \n ", weight_S_to_L4);
        % make sure weight doesn't cross limits
        if weight_S_to_L4 <= 0.001
            weight_S_to_L4 = 0.0001;
        end
        if weight_S_to_L4 >= 0.4
            weight_S_to_L4 = 0.4;
        end

        %-------- connections from D to L4
        voltage_L4_due_to_D = g_t4.*(weight_D_to_L4*spike_train_of_thalamus_D);
        [voltage_L4_due_to_D, spike_train_L4_due_to_D] =  decrease_voltage_for_20ms_after_spike(voltage_L4_due_to_D);
        
        arr_of_delta_ts = get_post_minus_pre_array(spike_train_of_thalamus_D, spike_train_L4_due_to_D); 
        for i=1:length(arr_of_delta_ts)
            if arr_of_delta_ts(1,i) > 0
                weight_D_to_L4 = weight_D_to_L4 * ( 1 + (a_LTP * exp(-arr_of_delta_ts(1,i)/tau_LTP)) );
            else
                weight_D_to_L4 = weight_D_to_L4 * ( 1 + (-a_LTD * exp(arr_of_delta_ts(1,i)/tau_LTD)) );
            end
        end
        % make sure weight doesn't cross limits
        if weight_D_to_L4 <= 0.001
            weight_D_to_L4 = 0.0001;
        end
        if weight_D_to_L4 >= 0.4
            weight_D_to_L4 = 0.4;
        end






        end % 1800 Stimulus - 9 m



        figure(11)
            plot(weights_from_S_to_L4);
            title("weights from S to L4");
        grid

        figure(12)
            plot(weights_from_D_to_L4);
            title("weights from D to L4");
        grid

        figure(13)
            plot(weights_from_SP_to_L4);
            title("weights from SP to L4");
        grid
end

function delta_ts_arr = get_post_minus_pre_array(presyn_arr, postsyn_arr)
    pre_syn_times = [];
    post_syn_times = [];
    for i=1:length(presyn_arr)
        if presyn_arr(1,i) == 1
            pre_syn_times = [pre_syn_times, i];
        end
    end

    for i=1:length(postsyn_arr)
        if postsyn_arr(1,i) == 1
            post_syn_times = [post_syn_times, i];
        end
    end

    smallest_length = min(length(pre_syn_times), length(post_syn_times));
    diff_arr = [];
    for i=1:smallest_length
        diff_arr = [diff_arr, post_syn_times(i) - pre_syn_times(i)];
    end

    delta_ts_arr = diff_arr;
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

function shifted_arr = shift_by_x(arr,x)
    % x is one, adjust with bad naming for now
    offset = 2;
    len = length(arr);
    new_arr = zeros(1, len);
    for i=offset:len
        new_arr(1,i) = arr(1, i-1); 
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