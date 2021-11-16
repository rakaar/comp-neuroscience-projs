function project
    stimulus = load('data_cn_project_iii_a17.mat', 'Stimulus').Stimulus;

    
 
    % Question - 1
    tau_plot = transpose(linspace(-50,50,101));
    autocorr_plot = [];
    for tau=-50:50
      autocorr = 0;
      n = 0;
      for i=1:20000
        if (i + tau >=1) & (i + tau <=20000)    
            autocorr = autocorr + (stimulus(1,i)*stimulus(1,i+tau));
            n = n+1;
        end
      end
      autocorr_plot = [autocorr_plot, autocorr/n];
    %   fprintf("for tau = %f, autocorr of stimulus is %f \n", tau, autocorr);
    end

    figure(101)
       plot(tau_plot, autocorr_plot);
    grid

     % Question-2
     all_spike_times = load('data_cn_project_iii_a17.mat', 'All_Spike_Times').All_Spike_Times;
   

     for neuron_num=1:4
         neuron_stimulus_rate = zeros(20000,1);
         for i=1:50
             spike_times = all_spike_times{neuron_num,i};
             len_of_spike_array = size(spike_times,2);
             for s=1:len_of_spike_array
                 % spike_times(s)*1000 as per this time, depends index of addition
                 index = fix(spike_times(s)*1000)+1;
                 neuron_stimulus_rate(index) = neuron_stimulus_rate(index) + 1;
 
             end
         end
 
         neuron_stimulus_rate = neuron_stimulus_rate./(50*0.001);
         time_axis = transpose(linspace(0,20,20000));
         
         disp(min(neuron_stimulus_rate));
         disp(max(neuron_stimulus_rate));
 
         figure(neuron_num)
             plot(time_axis, neuron_stimulus_rate);
         grid
 
     end


    % Question - 3
    colors = ['r', 'b', 'g', 'c'];
    bin_sizes = [0.01,0.02,0.05, 0.1, 0.2, 0.5];
    all_spike_times = load('data_cn_project_iii_a17.mat', 'All_Spike_Times').All_Spike_Times;
    
        for b=1:6
            mean_plot = zeros(4, 20*(1/bin_sizes(b)));
            variance_plot = zeros(4, 20*(1/bin_sizes(b)));
    
            for neuron_num=1:4
                sample = zeros(50, 20*(1/bin_sizes(b)));
                for i=1:50
                    spike_times = all_spike_times{neuron_num,i};
                    len_of_spike_array = size(spike_times,2);
                    
                    for s=1:len_of_spike_array
                        index = fix(spike_times(s)*(1/bin_sizes(b)))+1;
                        sample(i,index) = sample(i,index) + 1;
        
                    end
        
                
        
                end
        
                for i=1:20*(1/bin_sizes(b))
                    mean_plot(neuron_num, i) = mean(sample(:,i));
                    variance_plot(neuron_num, i) = var(sample(:,i));
                end
            end
            
            figure(30+b)
                 hold on
                    for k=1:4
                            scatter(mean_plot(k,:), variance_plot(k,:),[], colors(k));
                    end
                    max_var = max(variance_plot, [], 'all');
                    yisx = linspace(0, max_var);
                    plot(yisx, yisx, 'color', 'k');
                hold off
            grid
        
        end
   
     
    
  
    
end
    