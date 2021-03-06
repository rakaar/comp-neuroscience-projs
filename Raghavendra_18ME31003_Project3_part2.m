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
  
    
        
        figure(430)
        hold on
            for i=1:4
                plot(linspace(1,100), spike_triggered_avg(:, i));
            end
        hold off
        grid
    
    Rxx = autocorr(Stimulus, 99);
    Css = zeros(100,100);
    Css = toeplitz(Rxx);
    for i=1:4
        h_t(100:-1:1,i) = (Css\spike_triggered_avg(:,i));
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

        [fitresult1, gof] = fit( x1, y1, ft, opts );
        plot(fitresult1,x1, y1);
        hold off
    grid


    figure(52)
        hold on
        [x2, y2] = prepareCurveData( x2, y2 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];

       
        [fitresult2, gof] = fit( x2, y2, ft, opts );
        plot(fitresult2,x2, y2);
        hold off
    grid


    figure(53)
        hold on
        [x3, y3] = prepareCurveData( x3, y3 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];

       
        [fitresult3, gof] = fit( x3, y3, ft, opts );
        plot(fitresult3,x3, y3);
        hold off
    grid


    figure(54)
        hold on
        [x4, y4] = prepareCurveData( x4, y4 );

        ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.5 0.5, 0.5];


        [fitresult4, gof] = fit( x4, y4, ft, opts );
        plot(fitresult4,x4, y4);
        hold off
    grid


    %% Question 6
    % predicton - y , f
    % y


   
    y1 = conv(stimulus(15001:20000), h_t(:,1));
    y1 = y1(1,1:5000);
    
    y2 = conv(stimulus(15001:20000), h_t(:,2));
    y2 = y2(1,1:5000);
    
    y3 = conv(stimulus(15001:20000), h_t(:,3));
    y3 = y3(1,1:5000);
    
    y4 = conv(stimulus(15001:20000), h_t(:,4));
    y4 = y4(1,1:5000);

    disp("????")
    disp(fitresult1);
    % disp(fitresult1.a); % .a, .b, .c

    % prediction - lambda
    pred1 = (fitresult1.a)./(1+exp(-(fitresult1.b).*(y1-(fitresult1.c))));
    figure(61)
        scatter(PSTH(1, 15001:20000),pred1);
    grid

    c1 = corrcoef(PSTH(1, 15001:20000), pred1);
    disp("R square for Neuron 1")
    disp(c1(2)^2)


    
    pred2 = (fitresult2.a)./(1+exp(-(fitresult2.b).*(y2-(fitresult2.c))));
    figure(62)
        scatter(PSTH(2, 15001:20000),pred2);
    grid
    
    c2 = corrcoef(PSTH(2, 15001:20000), pred2);
    disp("R square for Neuron 2")
    disp(c2(2)^2)

    pred3 = (fitresult3.a)./(1+exp(-(fitresult3.b).*(y3-(fitresult3.c))));
    figure(63)
        scatter(PSTH(3, 15001:20000),pred3);
    grid

    c3 = corrcoef(PSTH(3, 15001:20000), pred3);
    disp("R square for Neuron 3")
    disp(c3(2)^2)
    
    pred4 = (fitresult4.a)./(1+exp(-(fitresult4.b).*(y4-(fitresult4.c))));
    figure(64)
        scatter(PSTH(4, 15001:20000),pred4);
    grid
    c4 = corrcoef(PSTH(4, 15001:20000), pred4);
    disp("R square for Neuron 4")
    disp(c4(2)^2)

    %% Question - 7
    % q = 0.1;
    q = [0 0.001 0.01 0.1 1 10 100];
    MI = zeros(4,5,length(q));

    for iter=1:5
        disp("i is ")
        disp(iter)
        disp("-------------------------------------------------")
        t = randperm(19901,8);

        for n = 1:4
            disp("n is ")
            disp(n)
            spike_segments = cell(8,50);
            for rep = 1:50
                spikes = All_Spike_Times{n,rep};
                for i = 1:8
                    spike_segments{i,rep} = spikes(spikes>=t(i)/1000&spikes<(t(i)+100)/1000);
                end
            end

            for m=1:length(q)
                disp('q is')
                disp(q(m))
                confusion_matrix = zeros(8,8);
                
                for i=1:8
                    for realiz=1:50
                        closest_stimuli_index = find_closest_stimuli_index(i, realiz, q, spike_segments);
                        % append to that stimulus
                        confusion_matrix(i, closest_stimuli_index) = confusion_matrix(i, closest_stimuli_index) + 1;
                    end
                end

                confusion_matrix = confusion_matrix/50;
                MI(n,iter,m) = MI(n,iter,m) + find_mutual_info(confusion_matrix);
                
            end

        end


    end

    

        
        for n = 1:4
            fprintf("neuron %f \f", n);
            avg_q_for_neuron = zeros(1, length(q));
                for i=1:length(q)
                    avg_q1 = sum(MI(n, :, i))/5;
                    avg_q_for_neuron(i) = avg_q1;
                    fprintf("for value of q - %f, MI are \n", q(i));
                    disp(MI(n,:,i));   
                end
                
                % figure(70 + n)                    
                %     plot(log10(q),avg_q_for_neuron);
                % grid
        end
        
    
end

function mutual_info = find_mutual_info(confusion_matrix)
    m = 0;
    for i=1:8
        for j=1:8
            if(confusion_matrix(i,j)~=0)
                m = m + ((confusion_matrix(i,j))* (1/8)) * log2((confusion_matrix(i,j))/sum(confusion_matrix(:,j)/8)) ; %p(y/x)* p(x) * log(p(x,y)/ p(x)p(y))  = p(y/x)* p(x) * log((p(y/x))/p(y)) [bcoz p(x,y)/p(x) = p(y/x)]  and p(y) =  for all x sum(p(x,y)) =  
            end
        end
    end

    mutual_info = m;
end

function closest_index = find_closest_stimuli_index(i, realiz, q, spike_segments)
    distances = zeros(1,8);
    for s=1:8
        for r = 1:50
            if (r == realiz && s == i)
                continue
            end
            distances(s) = distances(s) + VPSDM(spike_segments{i,realiz},spike_segments{s,r},q);
        end
    end

    [~,k] = min(distances);

    closest_index = k;
end

function d=VPSDM(tli,tlj,q)
    nspi=length(tli);
    nspj=length(tlj);
    
    if q==0
       d=abs(nspi-nspj);
       return
    elseif q==Inf
       d=nspi+nspj;
       return
    end
    
    scr=zeros(nspi+1,nspj+1);
    scr(:,1)=(0:nspi)';
    scr(1,:)=(0:nspj);
    if(nspi && nspj)
       for i=2:nspi+1
          for j=2:nspj+1
             scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+q*abs(tli(i-1)-tlj(j-1))]);
          end
       end
    end
    d=scr(nspi+1,nspj+1);
end