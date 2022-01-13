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
        voltage_L4 = [];

        % for S,D,SP4 neuron
        S_xe = [];
        S_xe = [S_xe, 0];
        S_xi = [];
        S_xi = [S_xi, 0];
        S_xr = [];
        S_xr = [S_xr, 1];


        D_xe = [];
        D_xe = [ D_xe, 0];
        D_xi = [];
        D_xi = [ D_xi, 0];
        D_xr = [];
        D_xr = [D_xr, 1];

        SP_xe = [];
        SP_xe = [SP_xe, 0];
        SP_xi = [];
        SP_xi = [SP_xi, 0];
        SP_xr = [];
        SP_xr = [SP_xr, 1];

        % tau values
        thalamus_tau_re = 0.9; thalamus_tau_ei=10; thalamus_tau_ir=5000;
        sp_tau_re = 0.9; sp_tau_ei = 27; sp_tau_ir = 5000;

        for i=1:1800 % 3 mins

            voltage_SP_for_single_stimulus = [];
            spike_for_S_for_single_stimulus = [];
            spike_for_D_for_single_stimulus = [];

            
            s_xe_trimmed_for_single_stimulus = [];
            d_xe_trimmed_for_single_stimulus = [];
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

         
            
            spike_for_S_for_single_stimulus = [spike_for_S_for_single_stimulus, spike_for_S];
            spike_for_D_for_single_stimulus = [spike_for_D_for_single_stimulus, spike_for_D];


            g_t_S1 = get_g_t(spike_for_S);
            g_t_D1 = get_g_t(spike_for_D);

            g_t_S1 = g_t_S1(1,1:50);
            g_t_D1 = g_t_D1(1,1:50);

            
            % synaptic resources for thalamic neurons
            for kkk=1:50
                Ms = 0; Md = 0;
                if spike_for_S(1,kkk) == 1
                    Ms = 1;
                end

                if spike_for_D(1,kkk) == 1
                    Md = 1;
                end
                current_S_xr = S_xr(1, length(S_xr));
                current_S_xe = S_xe(1, length(S_xe));
                current_S_xi = S_xi(1, length(S_xi));

                S_xr = [S_xr, current_S_xr + ((-Ms * (current_S_xr/thalamus_tau_re)) + (current_S_xi/thalamus_tau_ir)  )];
                S_xe = [S_xe, current_S_xe + ((Ms*(current_S_xr/thalamus_tau_re)) - (current_S_xe/thalamus_tau_ei))];
                S_xi = [S_xi, current_S_xi + (current_S_xe/thalamus_tau_ei) - (current_S_xi/thalamus_tau_ir)];
            
                current_D_xr = D_xr(1, length(D_xr));
                current_D_xe = D_xe(1, length(D_xe));
                current_D_xi = D_xi(1, length(D_xi));
            
                D_xr = [D_xr, current_D_xr + ((-Ms * (current_D_xr/thalamus_tau_re)) + (current_D_xi/thalamus_tau_ir)  )];
                D_xe = [D_xe, current_D_xe + ((Ms*(current_D_xr/thalamus_tau_re)) - (current_D_xe/thalamus_tau_ei))];
                D_xi = [D_xi, current_D_xi + (current_D_xe/thalamus_tau_ei) - (current_D_xi/thalamus_tau_ir)];
            end
            

            s_xe_trimmed = S_xe(1,length(S_xe)-49:length(S_xe));
            d_xe_trimmed = D_xe(1, length(D_xe)-49:length(D_xe));

            s_xe_trimmed_for_single_stimulus = [s_xe_trimmed_for_single_stimulus, s_xe_trimmed];
            d_xe_trimmed_for_single_stimulus = [d_xe_trimmed_for_single_stimulus, d_xe_trimmed];


            v_sp = weight_S_to_SP*shift_1(g_t_S1).*shift_1(spike_for_S).*s_xe_trimmed + weight_D_to_SP*shift_1(g_t_D1).*shift_1(spike_for_D).*d_xe_trimmed;
            voltage_SP = [voltage_SP, v_sp];
            voltage_SP_for_single_stimulus = [voltage_SP_for_single_stimulus, v_sp];


            % 250 ms gap
            spike_for_S = non_homo_poison(0.5, 1000);
            spike_for_D = non_homo_poison(0.5, 1000);

            spike_for_S = spike_for_S(1,1:250);
            spike_for_D = spike_for_D(1,1:250);

            spike_for_S_for_single_stimulus = [spike_for_S_for_single_stimulus, spike_for_S];
            spike_for_D_for_single_stimulus = [spike_for_D_for_single_stimulus, spike_for_D];

            % synaptic resources for thalamic neurons
            for kkk=1:250
                Ms = 0; Md = 0;
                if spike_for_S(1,kkk) == 1
                    Ms = 1;
                end

                if spike_for_D(1,kkk) == 1
                    Md = 1;
                end
                current_S_xr = S_xr(1, length(S_xr));
                current_S_xe = S_xe(1, length(S_xe));
                current_S_xi = S_xi(1, length(S_xi));

                S_xr = [S_xr, current_S_xr + ((-Ms * (current_S_xr/thalamus_tau_re)) + (current_S_xi/thalamus_tau_ir)  )];
                S_xe = [S_xe, current_S_xe + ((Ms*(current_S_xr/thalamus_tau_re)) - (current_S_xe/thalamus_tau_ei))];
                S_xi = [S_xi, current_S_xi + (current_S_xe/thalamus_tau_ei) - (current_S_xi/thalamus_tau_ir)];
            
                current_D_xr = D_xr(1, length(D_xr));
                current_D_xe = D_xe(1, length(D_xe));
                current_D_xi = D_xi(1, length(D_xi));
            
                D_xr = [D_xr, current_D_xr + ((-Ms * (current_D_xr/thalamus_tau_re)) + (current_D_xi/thalamus_tau_ir)  )];
                D_xe = [D_xe, current_D_xe + ((Ms*(current_D_xr/thalamus_tau_re)) - (current_D_xe/thalamus_tau_ei))];
                D_xi = [D_xi, current_D_xi + (current_D_xe/thalamus_tau_ei) - (current_D_xi/thalamus_tau_ir)];
            end

            g_t_S2 = get_g_t(spike_for_S);
            g_t_D2 = get_g_t(spike_for_D);


            g_t_S2 = g_t_S2(1,1:250);
            g_t_D2 = g_t_D2(1,1:250);
            
            s_xe_trimmed = S_xe(1,length(S_xe)-249:length(S_xe));
            d_xe_trimmed = D_xe(1, length(D_xe)-249:length(D_xe));

            s_xe_trimmed_for_single_stimulus = [s_xe_trimmed_for_single_stimulus, s_xe_trimmed];
            d_xe_trimmed_for_single_stimulus = [d_xe_trimmed_for_single_stimulus, d_xe_trimmed];


            
            v_sp = weight_S_to_SP*shift_1(g_t_S2).*shift_1(spike_for_S).*s_xe_trimmed + weight_D_to_SP*shift_1(g_t_D2).*shift_1(spike_for_D).*d_xe_trimmed;
            voltage_SP = [voltage_SP, v_sp];
            voltage_SP_for_single_stimulus = [voltage_SP_for_single_stimulus, v_sp];
            
            g_t_S3 = cat(2, g_t_S1, g_t_S2);
            g_t_D3 = cat(2, g_t_D1, g_t_D2);


          
            [not_needed_variable spike_train_SP_indiv_stimulus] = decrease_voltage_for_20ms_after_spike(voltage_SP_for_single_stimulus);
            g_sp = get_g_t(spike_train_SP_indiv_stimulus);
            g_sp = g_sp(1,1:300);
            
            % synaptic resources for SP
            for kkk=1:300
                Msp = 0;
                if spike_train_SP_indiv_stimulus(1,kkk) == 1;
                    Msp = 1;
                end

                current_SP_xr = SP_xr(1, length(SP_xr));
                current_SP_xe = SP_xe(1, length(SP_xe));
                current_SP_xi = SP_xi(1, length(SP_xi));

                SP_xr = [SP_xr, current_SP_xr + ((-Msp * (current_SP_xr/sp_tau_re)) + (current_SP_xi/sp_tau_ir)  )];
                SP_xe = [SP_xe, current_SP_xe + ((Msp*(current_SP_xr/sp_tau_re)) - (current_SP_xe/sp_tau_ei))];
                SP_xi = [SP_xi, current_SP_xi + (current_SP_xe/sp_tau_ei) - (current_SP_xi/sp_tau_ir)];

            end
            sp_xe_trimmed = SP_xe(1, length(SP_xe)-299:length(SP_xe));
            v_l4 = weight_S_to_L4* shift_1(g_t_S3).*shift_1(spike_for_S_for_single_stimulus).*s_xe_trimmed_for_single_stimulus + weight_D_to_L4*shift_1(g_t_D3).*shift_1(spike_for_D_for_single_stimulus).*d_xe_trimmed_for_single_stimulus + weight_SP_to_L4*shift_1(g_sp).*shift_1(spike_train_SP_indiv_stimulus).*sp_xe_trimmed;
            voltage_L4 = [voltage_L4, v_l4];
        end % end of for all stimulus

        [v_sp_modified, spike_train_for_SP] = decrease_voltage_for_20ms_after_spike(voltage_SP);
        for k=1:length(v_sp_modified)
            v_sp_modified(1,i) = v_sp_modified(1,i) * 0.9;
        end

        [v_l4_modified, spike_train_for_L4] = decrease_voltage_for_20ms_after_spike(voltage_L4);
        for k=1:length(v_l4_modified)
            v_l4_modified(1,i) = v_l4_modified(1,i) * 0.9;
        end

        figure(4)
            subplot(3,1,1)
            plot(S_xe);
            title('S xe')

            subplot(3,1,2)
            plot(D_xe);
            title('D xe')

            subplot(3,1,3)
            plot(SP_xe);
            title('SP xe')
        grid

        figure(3)
            subplot(2, 1, 1)
            plot(v_l4_modified)
            title('v L4')
            
            subplot(2, 1, 2)
            plot(spike_train_for_L4)
            title('spikes L4')
        grid


        figure(1)
            subplot(2, 1, 1)
            plot(v_sp_modified)
            title('v sp modified')
            
            subplot(2, 1, 2)
            plot(spike_train_for_SP)
            title('spike train SP')
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
        if voltage_values(1,i) < 0.03
            voltage_after_decreasing(1,i) = voltage_values(1,i);
            continue
        end

        if voltage_values(1,i) >= 0.03
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