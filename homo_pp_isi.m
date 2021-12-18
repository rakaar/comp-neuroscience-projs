function project
    homo_pp = adarsh_spk_gen_poss_nonhom(50, 0, 1, 0.001);
    figure(101)
        plot(homo_pp);
        title('homo');
    grid

    % isi
    isi = [];
    latest_spike_time = 0;
    inter_spike_value = 0;
    disp(size(homo_pp))

    for i=1:1000
        if homo_pp(1,i) == 0
            inter_spike_value = inter_spike_value + 1;
        else
            isi = [isi, inter_spike_value];
            inter_spike_value = 0;
            latest_spike_time = i;
        end    
    end
    isi(1) = [];
    figure(102)
        hist(isi);
        title('isi');
    grid


    % second way to test
    spike_times = [];
    for i=1:1000
        if homo_pp(1,i) == 1
            spike_times = [spike_times, i];
        end
    end
    isi = [];
    for i=1:length(spike_times)
        if i+1 <= length(spike_times)
            isi = [isi, spike_times(1,i+1) - spike_times(1,i)]; 
        end
    end
    figure(103)
        hist(isi);
        title('isi');
    grid

end

function spikeMat = adarsh_spk_gen_poss_nonhom(fr,tin,tout,dt)
    spk_mat=[];
    tSim=tout-tin;
    nBins=floor(tSim/dt);
    fprintf("nbins %d", nBins);
    spikeMat = zeros(1, nBins);
    rng('shuffle');
    for i=1:nBins
        if rand<=1-exp(-1*fr*dt);
            spikeMat(1,i)=1;
        else
            spikeMat(1,i)=0;
        end
    end
    % spk_mat=find(spikeMat==1);
    end