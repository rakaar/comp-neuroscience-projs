function project
    lambda_t = randi([40,50], 1,1000);
  testing_poison = non_homo_poison(lambda_t, 1000);

  figure(1)
    plot(testing_poison);
    title('poison process');
  grid

  % testing inter spike interval
  isi = [];
  latest_spike_time = 0;
  inter_spike_value = 0;
  for i=1:1000
    if testing_poison(1,i) == 0
        inter_spike_value = inter_spike_value + 1;
    else
        isi = [isi, inter_spike_value];
        inter_spike_value = 0;
        latest_spike_time = i;
    end    
  end

  figure(21)
    plot(isi);
    title('isi');
  grid

  % second way to testing isi
  spike_times = [];
  spike_times = [spike_times, 0]; 
  for i=1:1000
    if testing_poison(1,i) == 1
        spike_times = [spike_times, i];
    end
  end
  isi = [];
  for i=1:length(spike_times)
    if i+1 <= length(spike_times)
        isi = [isi, spike_times(1,i+1) - spike_times(1,i)]; 
    end
  end
  figure(22)
    plot(isi);
    title('isi');
  grid



  % testing adarsh code
  
    adarsh_non_homo_poisson = adarsh_spk_gen_poss_nonhom(lambda_t,0,1,0.001); 

    figure(41)
        plot(adarsh_non_homo_poisson);
    grid

    % isi
     % testing inter spike interval
  isi = [];
  latest_spike_time = 0;
  inter_spike_value = 0;
  for i=1:1000
    if adarsh_non_homo_poisson(1,i) == 0
        inter_spike_value = inter_spike_value + 1;
    else
        isi = [isi, inter_spike_value];
        inter_spike_value = 0;
        latest_spike_time = i;
    end    
  end

  figure(42)
    plot(isi);
    title('isi');
  grid

  % second way to testing isi
  spike_times = [];
  spike_times = [spike_times, 0]; 
  for i=1:1000
    if adarsh_non_homo_poisson(1,i) == 1
        spike_times = [spike_times, i];
    end
  end
  isi = [];
  for i=1:length(spike_times)
    if i+1 <= length(spike_times)
        isi = [isi, spike_times(1,i+1) - spike_times(1,i)]; 
    end
  end
  figure(43)
    plot(isi);
    title('isi');
  grid


end % end of proj


function spike_train = non_homo_poison(spike_rate, period)
    data = zeros(1, period);
    rng('shuffle');
    for i=1:period
            x = rand;
            if x < 1 - exp(-(spike_rate(1,i) * 0.001))
                data(1,i) = 1;
            end
    end

    spike_train = data;
end

function spikeMat = adarsh_spk_gen_poss_nonhom(fr,tin,tout,dt)
    spk_mat=[];
    tSim=tout-tin;
    nBins=floor(tSim/dt);
    fprintf("nbins %d", nBins);
    spikeMat = zeros(1, nBins);
    rng('shuffle');
    for i=1:nBins
        if rand<=1-exp(-1*fr(1,i)*dt);
            spikeMat(1,i)=1;
        else
            spikeMat(1,i)=0;
        end
    end
    % spk_mat=find(spikeMat==1);
    end