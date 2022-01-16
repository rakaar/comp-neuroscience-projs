function project 

    % vars
    weight_S_to_SP = 0.2; weight_D_to_SP = 0.2;
    weight_SP_to_L4 = 0.11; weight_S_to_L4 = 0.02; weight_D_to_L4 = 0.02;

    % all variables
    stimulus = []; % 1 for deviant, 0 for Standard
    voltage_S = []; spikes_S = [];
    voltage_D = []; spikes_D = [];
    voltage_SP = []; spikes_SP = [];

    % resources
    S_xr = []; S_xe = []; S_xi = [];
    D_xr = []; D_xe = []; D_xi = [];
    SP_xr = []; SP_xe = []; SP_xi = [];

    % initialize resources
    S_xr = [S_xr, 1]; D_xr = [D_xr, 1]; SP_xr = [SP_xr, 1];
    S_xe = [S_xe, 0]; D_xe = [D_xe, 0]; SP_xe = [SP_xe, 0];
    S_xi = [S_xi, 0]; D_xi = [D_xi, 0]; SP_xi = [SP_xi, 0]; 


    thalamus_tau_re = 0.9; thalamus_tau_ei = 10; thalamus_tau_ir = 5000;
    tau_SP_re = 0.9; tau_SP_ei = 27; tau_SP_ir = 5000;
    % voltage_L4 = []; spikes_L4 = [];
   
    NUMBER_OF_ITERS = 900;
    for i=1:NUMBER_OF_ITERS

        % first 50 ms stimulus - S or D
        stimulus_decider = randi([1,100]);
        if stimulus_decider >= 90
            % Deviant stimulus
            spikes_S_50ms = non_homo_poison(2.5, 50);
            spikes_D_50ms = non_homo_poison(10, 50);
            stimulus = [stimulus, 1];
        else 
            % Standard stimulus
            spikes_S_50ms = non_homo_poison(10, 50);
            spikes_D_50ms = non_homo_poison(2.5, 50);
            stimulus = [stimulus, 0];
        end

        spikes_S = [spikes_S, spikes_S_50ms];
        spikes_D = [spikes_D, spikes_D_50ms];

        % kernels
        g_t_S_50ms = get_g_t(spikes_S_50ms);
        g_t_S_50ms = g_t_S_50ms(1,1:50);

        g_t_D_50ms = get_g_t(spikes_D_50ms);
        g_t_D_50ms = g_t_D_50ms(1,1:50);

        % resources depletion for S and D
        for i=1:50
            Ms = 0; Md = 0;
            if spikes_S_50ms(1,i) == 1
                Ms = 1;
            end

            if spikes_D_50ms(1,i) == 1
                Md = 1;
            end
            % for S 
            current_S_xe = S_xe(1, length(S_xe));
            current_S_xr = S_xr(1, length(S_xr));
            current_S_xi = S_xi(1, length(S_xi));

            S_xr = [S_xr, update_xr(Ms, current_S_xr, current_S_xi, thalamus_tau_re, thalamus_tau_ir)];
            S_xe = [S_xe, update_xe(Ms, current_S_xr, current_S_xe, thalamus_tau_re, thalamus_tau_ei)];
            S_xi = [S_xi, update_xi(current_S_xe, current_S_xi, thalamus_tau_ei, thalamus_tau_ir)];

            % for D
            current_D_xe = D_xe(1, length(D_xe));
            current_D_xr = D_xr(1, length(D_xr));
            current_D_xi = D_xi(1, length(D_xi));

            D_xr = [D_xr, update_xr(Md, current_D_xr, current_D_xi, thalamus_tau_re, thalamus_tau_ir)];
            D_xe = [D_xe, update_xe(Md, current_D_xr, current_D_xe, thalamus_tau_re, thalamus_tau_ei)];
            D_xi = [D_xi, update_xi(current_D_xe, current_D_xi, thalamus_tau_ei, thalamus_tau_ir)];
        end
        S_xe_trimmed_50ms = S_xe(1, length(S_xe)-49:length(S_xe));
        D_xe_trimmed_50ms = D_xe(1, length(D_xe)-49:length(D_xe));

        % calculate voltage of SP 50ms
        voltage_SP_50ms = weight_S_to_SP*shift_1(g_t_S_50ms).*shift_1(spikes_S_50ms).*S_xe_trimmed_50ms  + weight_D_to_SP*shift_1(g_t_D_50ms).*shift_1(spikes_D_50ms).*D_xe_trimmed_50ms;
        voltage_SP = [voltage_SP, voltage_SP_50ms];


        % for 250 ms gap
        spikes_S_250ms = non_homo_poison(0.5, 250);
        spikes_D_250ms = non_homo_poison(0.5, 250);

        spikes_S = [spikes_S, spikes_S_250ms];
        spikes_D = [spikes_D, spikes_D_250ms];

        % kernels
        g_t_S_250ms = get_g_t(spikes_S_250ms);
        g_t_S_250ms = g_t_S_250ms(1,1:250);

        g_t_D_250ms = get_g_t(spikes_D_250ms);
        g_t_D_250ms = g_t_D_250ms(1,1:250);

        % resources depletion for S and D
        for i=1:250
            Ms = 0; Md = 0;
            if spikes_S_250ms(1,i) == 1
                Ms = 1;
            end

            if spikes_D_250ms(1,i) == 1
                Md = 1;
            end
            % for S 
            current_S_xe = S_xe(1, length(S_xe));
            current_S_xr = S_xr(1, length(S_xr));
            current_S_xi = S_xi(1, length(S_xi));

            S_xr = [S_xr, update_xr(Ms, current_S_xr, current_S_xi, thalamus_tau_re, thalamus_tau_ir)];
            S_xe = [S_xe, update_xe(Ms, current_S_xr, current_S_xe, thalamus_tau_re, thalamus_tau_ei)];
            S_xi = [S_xi, update_xi(current_S_xe, current_S_xi, thalamus_tau_ei, thalamus_tau_ir)];

            % for D
            current_D_xe = D_xe(1, length(D_xe));
            current_D_xr = D_xr(1, length(D_xr));
            current_D_xi = D_xi(1, length(D_xi));

            D_xr = [D_xr, update_xr(Md, current_D_xr, current_D_xi, thalamus_tau_re, thalamus_tau_ir)];
            D_xe = [D_xe, update_xe(Md, current_D_xr, current_D_xe, thalamus_tau_re, thalamus_tau_ei)];
            D_xi = [D_xi, update_xi(current_D_xe, current_D_xi, thalamus_tau_ei, thalamus_tau_ir)];
        end
        S_xe_trimmed_250ms = S_xe(1, length(S_xe)-249:length(S_xe));
        D_xe_trimmed_250ms = D_xe(1, length(D_xe)-249:length(D_xe));

        % calculate voltage of SP 250ms
        voltage_SP_250ms = weight_S_to_SP*shift_1(g_t_S_250ms).*shift_1(spikes_S_250ms).*S_xe_trimmed_250ms  + weight_D_to_SP*shift_1(g_t_D_250ms).*shift_1(spikes_D_250ms).*D_xe_trimmed_250ms;
        voltage_SP = [voltage_SP, voltage_SP_250ms];


    end % -------------- END OF STIMULUS --------


    [calculated_voltage_SP spikes_SP] = calculate_voltage_and_spikes(voltage_SP); 
    figure(1)
        subplot(3, 1, 1)
        plot(voltage_SP);
        title('raw voltage SP');

        subplot(3, 1, 2)
        plot(spikes_SP);
        title('spikes SP')

        subplot(3, 1, 3)
        plot(calculated_voltage_SP);
        title('calculated voltage SP');
    grid

    figure(2)
        subplot(3,1,1)
        stem(stimulus);
        title('stimuls 1 = deviant, 0 = standard')

        subplot(3,1,2)
        stem(spikes_D);
        title('spikes D');

        subplot(3,1,3)
        stem(spikes_S);
        title('spikes S');
    grid
end % END OF PROJECT ------------



function new_xr = update_xr(M, x_r, x_i, tau_re, tau_ir)
    new_xr = x_r + (-M*(x_r/tau_re)) + (x_i/tau_ir);
end

function new_xe = update_xe(M, x_r, x_e, tau_re, tau_ei)
    new_xe = x_e + (M*(x_r/tau_re)) - (x_e/tau_ei);
end

function new_xi = update_xi(x_e, x_i, tau_ei, tau_ir)
    new_xi = x_i + (x_e/tau_ei) - (x_i/tau_ir);
end

function [calculated_voltage spikes] = calculate_voltage_and_spikes(original_voltage)
    calculated_voltage = original_voltage;
    spikes = zeros(1, length(original_voltage));
    threshold = 0.005;
    t_spike = 0; tau_ref = 2;beta = 5;
    % for i=1:length(calculated_voltage)
    %     if calculated_voltage(1,i) >= threshold
    %         spikes(1,i) = 1;
    %         t_spike = i;
    %         for j=i:i+19
    %             if j <= length(calculated_voltage)
    %                 calculated_voltage(1,i) = calculated_voltage(1,i)*(1-(beta*exp(-(j-t_spike)/tau_ref)));
    %             end    
    %         end
    %     end 

    %     i = j; % start where the decreasing ends
    % end

    i=1;
    while i <= length(calculated_voltage)
        if calculated_voltage(1,i) >= threshold
            spikes(1,i) = 1;
            disp(i);
            t_spike = i;
            for j=i:i+19
                if j <= length(calculated_voltage)
                    calculated_voltage(1,i) = calculated_voltage(1,i)*(1-(beta*exp(-(j-t_spike)/tau_ref)));
                end    
            end
            i = j; % start where the decreasing ends
        else
            i = i + 1;
        end 
    end

    % decreasing 10% 
    for i=1:length(calculated_voltage)
        calculated_voltage(1,i) = calculated_voltage(1,i) * 0.9;
    end
end

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
