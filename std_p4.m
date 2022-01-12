function project

        % weights from S,D to SP    
        weight_S_to_SP = 0.2;
        weight_D_to_SP = 0.2;

        weight_SP_to_L4 = 0.11;
        weight_S_to_L4 = 0.02;
        weight_D_to_L4 = 0.02;

        a_LTP = 0.015;
        tau_LTP = 13;

        a_LTD = 0.021;
        tau_LTD = 20;
        
        voltage_SP = [];

        for i=1:600 % 3 mins
            stimulus_decider = randi([1,100]);
            spike_train_of_thalamus_S = [];
            spike_train_of_thalamus_D = [];


            if stimulus_decider >= 90
                % Deviant stimulus
                spike_for_S = non_homo_poison(2.5, 1000);
                spike_for_D = non_homo_poison(10, 1000);
            else 
                % Standard stimulus
                spike_for_S = non_homo_poison(10, 1000);
                spike_for_D = non_homo_poison(2.5, 1000);
            end

            % 50 ms stimulus
            spike_for_S = spike_for_S(1,1:50);
            spike_for_D = spike_for_D(1,1:50);

         



            g_t_S = get_g_t(spike_for_S);
            g_t_D = get_g_t(spike_for_D);

            g_t_S = g_t_S(1,1:50);
            g_t_D = g_t_D(1,1:50);

            v_sp = weight_S_to_SP*shift_1(g_t_S).*shift_1(spike_for_S) + weight_D_to_SP*shift_1(g_t_D).*shift_1(spike_for_D);
            voltage_SP = [voltage_SP, v_sp];

            % 250 ms gap
            spike_for_S = non_homo_poison(0.5, 1000);
            spike_for_D = non_homo_poison(0.5, 1000);

            spike_for_S = spike_for_S(1,1:250);
            spike_for_D = spike_for_D(1,1:250);


            g_t_S = get_g_t(spike_for_S);
            g_t_D = get_g_t(spike_for_D);


            g_t_S = g_t_S(1,1:250);
            g_t_D = g_t_D(1,1:250);
            


          v_sp = weight_S_to_SP*shift_1(g_t_S).*shift_1(spike_for_S) + weight_D_to_SP*shift_1(g_t_D).*shift_1(spike_for_D);
          voltage_SP = [voltage_SP, v_sp];
            

        end % end of for 1:1800

        [v_sp_modified, spike_train_for_SP] = decrease_voltage_for_20ms_after_spike(voltage_SP);
        for k=1:length(v_sp_modified)
            v_sp_modified(1,i) = v_sp_modified(1,i) * 0.9;
        end

        figure(1)
            plot(v_sp_modified);
        grid

        figure(2)
            plot(spike_train_for_SP);
        grid

end % end of project

function spike_train = non_homo_poison(spike_rate, period)
    data = zeros(1, period);
    rng('shuffle');
    for i=1:period
            x = rand;
            if x <= 1 - exp(-(spike_rate * 0.001))
                data(1,i) = 1;
            end
    end

    spike_train = data;
end

function k_t = get_k_t(spike_train)
    t_spike = 0;
    k_t = zeros(1,length(spike_train));
    for i=1:length(spike_train)
        if spike_train(1,i) == 1
            t_spike = i;
        end
        k_t(1,i) = exp(-(i-t_spike)/10);
    end    
end

function g_t = get_g_t(spike_train)
    kernel_kt = shift_1(get_k_t(spike_train));
    g_t = conv(kernel_kt, spike_train);
end

function shifted_arr = shift_1(arr)
    len = length(arr);
    new_arr = zeros(1, len);
    for i=2:len
        new_arr(1,i) = arr(1, i-1); 
    end

    shifted_arr = new_arr;
end

function  [new_voltage_values spike_train] = decrease_voltage_for_20ms_after_spike(voltage_values)
    voltage_after_decreasing = zeros(1, length(voltage_values));
    nearest_spike_time = 0;

    actual_spike = zeros(1, length(voltage_values));

    % variables for v_delta
    beta = 5; tau = 2;

    for i=1:length(voltage_values)
        if voltage_values(1,i) < 0.05
            voltage_after_decreasing(1,i) = voltage_values(1,i);
            continue
        end

        if voltage_values(1,i) >= 0.05
            nearest_spike_time = i;
            actual_spike(1,i) = 1;

            for j=i:i+19
                if j < length(voltage_values)
                    voltage_after_decreasing(1,j) = voltage_values(1,j)*(1 - ( beta * exp(-(j-nearest_spike_time)/tau) ))  ; 
                end
            end
        end

        i = j;
        
    end


    new_voltage_values = voltage_after_decreasing;
    spike_train = actual_spike;
end