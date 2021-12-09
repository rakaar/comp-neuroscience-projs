function project

    %% non homo
    
    % a poisson distribution with just 0s and 1s that is non homogeneous = binomial distrubutiion with varying probability
    random_small_prob = randi([40,50], 100,1000)/1000;
    nonhomopp = binornd(1, random_small_prob, 100,1000);

    disp(max(nonhomopp, [], 'all'))
    
    % generate a homogenous poison process with different repitions
    % reps = [10, 20,40,80]
    % TODO - RMSE of exactly what?? PSTH is 1x1000, the reps  are 100x1000, how rmse?
    
    %% going with homog poison
    random_nums = randi([1, 1000], 100,1000);
    poisson_spike = zeros(100,1000);
    for j=1:100
        for i=1:1000
            if random_nums(j,i) <= 45
                poisson_spike(j,i) = 1;
            end 
        end
    end
    
    PSTH = zeros(1,1000);
    for j=1:1000

        for i=1:100
            PSTH(1,j) = PSTH(1,j) + poisson_spike(i,j); 
        end
        PSTH(1,j) = PSTH(1,j)/(100*0.001);

    end

    plot(linspace(1,1000,1000), PSTH);
end

function count = get_occurance(spike_train, ele)
    index = 0;
    for i=1:length(spike_train)
        if spike_train(1,i) == ele
            index = index + 1;
        end
    end
    count = index;
    
end
% i don't understand non homogenous pp for now. going with homogenous
% function y = nonhomopp(intens,T)
%     % example of generating a 
%     % nonhomogeneousl poisson process on [0,T] with intensity function intens
    
%     x = 0:.1:T;
%     m = intens(x);
%     m2 = max(m); % generate homogeneouos poisson process
%     u = rand(1,ceil(1.5*T*m2));
%     y = cumsum(-(1/m2)*log(u)); %points of homogeneous pp
%     y = y(y<T); n=length(y); % select those points less than T
%     m = intens(y); % evaluates intensity function
%     y = y(rand(1,n)<m/m2); % filter out some points
%     histogram(y,10)
% end    