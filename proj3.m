function project
    stimulus = load('data_cn_project_iii_a17.mat', 'Stimulus').Stimulus;

    % Question - 1
    for tau=-50:50
      autocorr = 0;
      for i=1:20000
        if (i + tau >=1) & (i + tau <=20000)    
            autocorr = autocorr + (stimulus(1,i)*stimulus(1,i+tau));
            autocorr = autocorr/20000; % m value = 20000

        end
      end

      fprintf("for tau = %f, autocorr of stimulus is %f \n", tau, autocorr);
    end
end