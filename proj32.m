function project
    
    %% Question-4
    load('data_cn_project_iii_a17.mat');


    spike_triggered_avg = zeros(4,100);
    stimulus = load('data_cn_project_iii_a17.mat', 'Stimulus').Stimulus;
    all_spike_times = load('data_cn_project_iii_a17.mat', 'All_Spike_Times').All_Spike_Times;
    
    for neuron_num=1:4
    % neuron_num = 1;
        total_num_of_spikes = 0;
        for i=1:50
            spike_times = all_spike_times{neuron_num,i};
            len_of_spike_array = size(spike_times,2);
            
            for j=1:len_of_spike_array
                if spike_times(j) <= 15
                    total_num_of_spikes = total_num_of_spikes + 1;
                    time_of_spike_in_ms = spike_times(j)*1000;

                    % disp("errr belo")
                    % disp(size(spike_triggered_avg));

                    for k=1:100
                    
                        if time_of_spike_in_ms-(101-k) >=1
                            spike_triggered_avg(neuron_num, k) = spike_triggered_avg(neuron_num, k) +  stimulus(1,floor(time_of_spike_in_ms-(101-k)));
                        end
                    end
                end
            end
        end
        spike_triggered_avg(neuron_num,:) = spike_triggered_avg(neuron_num,:)./total_num_of_spikes; %  4 x 100
    end
   
    spike_triggered_avg = transpose(spike_triggered_avg); % 100 x 4
  
    for i=1:4
        figure(400+i)
            plot(linspace(1,100), spike_triggered_avg(:, i));
        grid
    end
    Rxx = autocorr(Stimulus, 99);
    Css = zeros(100,100);
    Css = toeplitz(Rxx);
    for i=1:4
        % h(i,100:-1:1)=(Css\sta(i,:)')';
        h_t(100:-1:1,i) = (Css\spike_triggered_avg(:,i));
        %  disp(toeplitz(autocorr(stimulus, 99)))
        % disp(spike_triggered_avg(:,i));
        figure(40+i)
            plot(linspace(1, 100), h_t(:,i));
        grid
        
        
    end

  
    %% Question - 5
    for i=1:4
        y_t(i,:) = conv(Stimulus(1, 1:15000), h_t(:, i));
        
    end
    disp(size(y_t)) % 4 x 15099

    PSTH = zeros(4, 20000);
    for i=1:4
        for j=1:50
            PSTH(i,:) = PSTH(i,:) + histcounts(All_Spike_Times{i,j}*1000,0:20000)*1000/50;
        end
    end
   
    for i=1:4
        y_t1(i, 1:15000) = y_t(i,1:15000);
    end
    % average the results in bins
    bin_size = 30;
    for i = 1:ceil(15000/bin_size)
        e = i*bin_size;
        if e>15000
            e = 15000;
        end
        x1(i) = mean( y_t(1, (1+(i-1)*bin_size):e) );
        x2(i) = mean( y_t(2, (1+(i-1)*bin_size):e) );
        x3(i) = mean( y_t(3, (1+(i-1)*bin_size):e) );
        x4(i) = mean( y_t(4, (1+(i-1)*bin_size):e) );
    
        y1(i) = mean( PSTH(1, (1+(i-1)*bin_size):e) );
        y2(i) = mean( PSTH(2, (1+(i-1)*bin_size):e) );
        y3(i) = mean( PSTH(3, (1+(i-1)*bin_size):e) );
        y4(i) = mean( PSTH(4, (1+(i-1)*bin_size):e) );
    end
    
    % xi = h(t) - 1 x 15000, yi = PSTHi - 1 x 15000
    

    figure(51)
        hold on
        [x1, y1] = prepareCurveData( x1, y1 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];

        [fitresult, gof] = fit( x1, y1, ft, opts );
        plot(fitresult,x1, y1);
        hold off
    grid


    figure(52)
        hold on
        [x2, y2] = prepareCurveData( x2, y2 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];

       
        [fitresult, gof] = fit( x2, y2, ft, opts );
        plot(fitresult,x2, y2);
        hold off
    grid


    figure(53)
        hold on
        [x3, y3] = prepareCurveData( x3, y3 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];

       
        [fitresult, gof] = fit( x3, y3, ft, opts );
        plot(fitresult,x3, y3);
        hold off
    grid


    figure(54)
        hold on
        [x4, y4] = prepareCurveData( x4, y4 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];


        [fitresult, gof] = fit( x4, y4, ft, opts );
        plot(fitresult,x4, y4);
        hold off
    grid



   

    
    
    



end