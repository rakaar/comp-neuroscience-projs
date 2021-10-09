

function project
    v = linspace(-80, 80);

    % defining first set of MLE variables
    g_ca = 4.4;
    g_k = 8;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.02;
    v1 = -1.2;
    v2 = 18;
    v3 = 2;
    v4 = 30;
    v5 = 2;
    v6 = 30;
    c = 20;
    
    
    % plotting null-clines,finding equilibrium points and their nature, quiver plots
    % MLE with first set of variables, i external = 0
    figure(1)
    hold on
        % w null-cline dw/dt = 0, w = f(v)
        w_null_cline = 0.5 * ( 1 + tanh((v-v3)/(v4)) );
        
        % v null cline, dv/dt = 0, w = f(v)
        m_infinity_v = 0.5 * ( 1 + tanh((v-v1)/(v2)) ); 
        denominator_v_null_cline = g_k * (v - v_k);
        numerator_v_null_cline = -g_ca * ( m_infinity_v.* (v-v_ca) ) - g_l * (v-v_l);
        
        v_null_cline = numerator_v_null_cline./denominator_v_null_cline;

        plot(v, 100*w_null_cline)
        plot(v, 100*v_null_cline)

      
        % point of intersection- method  1 - iterating in vectors
        for cnt=1:numel(v)
            s1 = v_null_cline(cnt);
            s2 = w_null_cline(cnt);
            zz = abs(s1-s2);
            if zz <= 0.01
               fprintf("v_eq,w_eq %.3f, %.3f \n", v(cnt), s1)
               v_eq = v(cnt)
               w_eq = s1
            end
        end
        
        % findng point of intersection - method 2 - using sym package of matlab(NR Method)
        F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)));  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
        starting_pt = [-60; 0.01];
        options = optimoptions('fsolve','Display','iter');
        [x,fval] = fsolve(F,starting_pt,options)
        disp("Equilibrium point for MLE with first set of variables, i external = 0")
        disp(x)
        


        % quiver plots
        [v_quiver,w_quiver] = meshgrid(linspace(-80,80,30), linspace(0,1,30));
        m_infinity_v_quiver = 0.5 * ( 1 + tanh((v_quiver-v1)/(v2)) ); 

        tau_w = 1./cosh((v_quiver-v3)/(2*v4));
        dv_dt = (1/c)*((-g_ca * ( m_infinity_v_quiver.*(v_quiver-v_ca) )) + (-g_k * ( w_quiver.*(v_quiver-v_k) )) + (-g_l * (v_quiver - v_l)));
        dw_dt = phi * (0.5 * ( 1 + tanh((v_quiver-v3)/(v4)) ) - w_quiver)./tau_w;
        quiver(v_quiver,100*w_quiver, dv_dt, 100*dw_dt, 10, 'color',[0 0 0]); % arrow length scaled 2 times for visibility
        
        
        % defining dv/dt and dw/dt again because this time we want expressions not vectors for the jacobian
        % dv/dt and dw_dt are functions of (v_var, w_var)
        syms v_var w_var
        dv_dt2 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var-v1)/(v2)) ))*(v_var-v_ca) )) + (-g_k * ( w_var*(v_var-v_k) )) + (-g_l * (v_var - v_l)));
        dw_dt2 = phi * (0.5 * ( 1 + tanh((v_var-v3)/(v4)) ) - w_var)/(1/cosh((v_var-v3)/(2*v4)));
        
        df1_dv = diff(dv_dt2, v_var);
        df1_dw = diff(dv_dt2, w_var);
        df2_dv = diff(dw_dt2, v_var);
        df2_dw = diff(dw_dt2, w_var);

        % jacobian matrix and their eigen values
        disp("jacobian matrix of MLE Iext=0")
        jacobian = [subs(df1_dv,{v_var,w_var},{v_eq, w_eq}) subs(df1_dw,{v_var,w_var},{v_eq, w_eq}); subs(df2_dv,{v_var,w_var},{v_eq, w_eq}) subs(df2_dw,{v_var,w_var},{v_eq, w_eq})  ];
        disp(double(jacobian));
        disp("Eigen values of MLE I ext=0, should be stable - both eigen values negative")
        eigen_values = double(eig(jacobian)) % we see that eigen values are negative, implying that equilibrium point is a stable point
        
        % running the equation from Equilibrium point, it stays there
        [t,r] = ode15s(@mle_diff_eqn,[0 100],[x(1) x(2)])
        plot(r(:,1),100*r(:,2));

    hold off
    grid

    
    
    
    
    % plotting for different values of phi
    % MLE Different set of variables for 
    figure(2)
        hold on
        % [t1,r1] = ode15s(@mle_diff_eqn,[0 10000],[x(1)+10 x(2)]) % sub threshold, phi = 0.02
        % plot4 = plot(r1(:,1),r1(:,2))

        
        + 30
        % [t1,r1] = ode15s(@mle_diff_eqn2,[0 10000],[x(1)+10 x(2)]) % sub threshold, phi = 0.04
        % plot5 = plot(r1(:,1),r1(:,2))

        

        % [t1,r1] = ode15s(@mle_diff_eqn3,[0 10000],[x(1)+10 x(2)]) % sub threshold, phi = 0.01
        % plot6 = plot(r1(:,1),r1(:,2))

        [t,r] = ode15s(@mle_diff_eqn3,[0 300],[x(1)+100 x(2)]); % supra threshold, phi = 0.01
        plot1 = plot(r(:,1), r(:,2))

        [t,r] = ode15s(@mle_diff_eqn,[0 300],[x(1)+100 x(2)]); % supra threshold, phi = 0.02
        plot2 = plot(r(:,1), r(:,2))

        [t,r] = ode15s(@mle_diff_eqn2,[0 300],[x(1)+100 x(2)]); % supra threshold, phi = 0.04
        plot3 = plot(r(:,1), r(:,2));

        legend([plot1, plot2, plot3], ["AP phi=0.01", "AP phi=0.02", "AP phi = 0.04"]);
       
        hold off
    grid

    % plotting for different values of V_initial to see sub and supra thresholds
    figure(3)
        hold on
            plots_depoloraization = [];
            labels = [];

            for step_size = 40:1:50
                [t,r] = ode15s(@mle_diff_eqn,[0 300],[x(1)+step_size x(2)]);
                plots_depoloraization = [plots_depoloraization, plot(r(:,1), r(:,2))];
                labels = [labels, strcat("V initial=",string(x(1) + step_size))];
            end
            legend(plots_depoloraization, labels)
        hold off
    grid

    % V max vs V initial
    % The plot shows the presence of threshold voltage, intial value above threshold voltage results in action potential kind behaviour
    figure(4)
        hold on
            v_initial = [];
            v_max = [];
            for step_size = 45:0.01:46
                [t,r] = ode15s(@mle_diff_eqn,[0 300],[x(1)+step_size x(2)]);
                v_initial = [v_initial, x(1)+step_size];
                v_max = [v_max, max(r(:,1))];
            end
        hold off

        plot(v_initial, v_max)
    grid

    % with steady external currents
    figure(5)
        hold on
            [t,r] = ode15s(@mle_diff_eqn,[0 300],[x(1) x(2)]); % equilibrium point, iext = 0
            plot1 = plot(r(:,1), r(:,2)) %nothing happens as we are at the equilibrium point

            [t,r] = ode15s(@mle_diff_eqn_with_i_ext_steady,[0 300],[x(1) x(2)]); % equilibrium point, iext = 86
            plot2 = plot(r(:,1), r(:,2)) % we get a limit cycle, we are starting from a point far from new equilibrium point, so can't predict via linearisation, derivates play the role

            [t,r] = ode15s(@mle_diff_eqn_with_i_ext_steady,[0 500],[-27.9, 0.17]); % non-equilibirum point, iext = 86
            plot3 = plot(r(:,1), r(:,2)) % we get a inward spiral, we are starting a point near new equilbirium point. And jacobian eigen values given complex with real part -ve, which confirms inward spiral

            legend([plot1, plot2, plot3], ["equilibrium pt, iext = 0", "equilibrium point, iext = 86", "non-equilibirum point, iext = 86"]);
        hold off

        
        % Calculations verifying the plots in figure(5)
        % finding - new equilibrium point when i ext = 86, its jacobian, its eigen values
        F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)) + 86);  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
        starting_pt = [-27; 0.10];
        options = optimoptions('fsolve','Display','iter');
        [x,fval] = fsolve(F,starting_pt,options)
        disp("Equilibrium pt with iext = 86")
        disp(x) % -27.9524, 0.1195
        v_eq3 = x(1)
        w_eq3 = x(2)

        % finding jacobian values
        syms v_var3 w_var3
        dv_dt3 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var3-v1)/(v2)) ))*(v_var3-v_ca) )) + (-g_k * ( w_var3*(v_var3-v_k) )) + (-g_l * (v_var3 - v_l)) + 86);
        dw_dt3 = phi * (0.5 * ( 1 + tanh((v_var3-v3)/(v4)) ) - w_var3)/(1/cosh((v_var3-v3)/(2*v4)));
        
        df1_dv3 = diff(dv_dt3, v_var3);
        df1_dw3 = diff(dv_dt3, w_var3);
        df2_dv3 = diff(dw_dt3, v_var3);
        df2_dw3 = diff(dw_dt3, w_var3);

        % jacobian matrix and their eigen values
        jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
        disp("eigen values for i ext = 86")
        eigen_values3 = double(eig(jacobian3)) % we see that eigen values are negative, implying that equilibrium point is a stable point
        % eigen values are complex, with real part negative => inward spiral
        % we get that inward spiral when we start near the equilibrium pt -27.9,0.17

    grid

    
    
    % i = 86
    % run backwards in time, find UPO,  null clines
    figure(6)
        hold on
             % plot for null clines with iext = 86
            % w null-cline dw/dt = 0, w = f(v)
            w_null_cline = 0.5 * ( 1 + tanh((v-v3)/(v4)) );
                
            % v null cline, dv/dt = 0, w = f(v)
            m_infinity_v = 0.5 * ( 1 + tanh((v-v1)/(v2)) ); 
            denominator_v_null_cline = g_k * (v - v_k);
            numerator_v_null_cline = -g_ca * ( m_infinity_v.* (v-v_ca) ) - g_l * (v-v_l) + 86;
            
            v_null_cline = numerator_v_null_cline./denominator_v_null_cline;

            
            plot(v, 100*w_null_cline)
            plot(v, 100*v_null_cline)
            
            
            [t,r] = ode15s(@mle_diff_eqn_with_i_ext_steady_backward_time,[0 300],[-27.9 0.17]);
            plot(r(:,1), 100*r(:,2)) 



            [t,r] = ode15s(@mle_diff_eqn_with_i_ext_steady,[0 300],[-37.98 12.60]);
            plot(r(:,1), 100*r(:,2))

            [t,r] = ode15s(@mle_diff_eqn_with_i_ext_steady,[0 300],[-30 0.15]);
            plot(r(:,1), 100*r(:,2))

            xlim([-60 0])
            ylim([0 50])

        hold off
    grid


    % other all i currents comparison
    iext = [80, 86, 90]
    for i=1:3
        disp("for i ext is ")
        F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)) + iext(i));  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
        starting_pt = [-27; 0.10];
        options = optimoptions('fsolve','Display','iter');
        [x,fval] = fsolve(F,starting_pt,options)
        disp(x) % -27.9524, 0.1195
        v_eq3 = x(1)
        w_eq3 = x(2)

        % finding jacobian values
        syms v_var3 w_var3
        dv_dt3 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var3-v1)/(v2)) ))*(v_var3-v_ca) )) + (-g_k * ( w_var3*(v_var3-v_k) )) + (-g_l * (v_var3 - v_l)) + iext(i));
        dw_dt3 = phi * (0.5 * ( 1 + tanh((v_var3-v3)/(v4)) ) - w_var3)/(1/cosh((v_var3-v3)/(2*v4)));
        
        df1_dv3 = diff(dv_dt3, v_var3);
        df1_dw3 = diff(dv_dt3, w_var3);
        df2_dv3 = diff(dw_dt3, v_var3);
        df2_dw3 = diff(dw_dt3, w_var3);

        % jacobian matrix and their eigen values
        jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
        eigen_values3 = double(eig(jacobian3)) 
    end

    % frequency of action potential vs current
    figure(10)
        rates_of_ap = [];
        i_ext = [];
        for i=80:100
            i_ext = [i_ext, i];
            [t r] = mle_solution_i_ext_set1(i);
            frequency_of_ap = 1/calculate_ap_time(r,t);
            rates_of_ap = [rates_of_ap, frequency_of_ap];
        end

        plot(i_ext, rates_of_ap);
    grid
    
    
    % To see what happens for equilibirum points as applied current changes
    % calculate equilibrium points from 80 to 100, plot them on (v,w) curve

    % -------------------------------------------------- Different set of MLE variables ------------------------------------------
    g_ca = 4;
    g_k = 8.0;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.0667;
    v1 = -1.2;
    v2 = 18;
    v3 = 12;
    v4 = 17.4;
    v5 = 12;
    v6 = 17.4;
    c = 20;

    % w null-cline dw/dt = 0, w = f(v)
    w_null_cline = 0.5 * ( 1 + tanh((v-v3)/(v4)) );
        
    % v null cline, dv/dt = 0, w = f(v)
    m_infinity_v = 0.5 * ( 1 + tanh((v-v1)/(v2)) ); 
    denominator_v_null_cline = g_k * (v - v_k);
    numerator_v_null_cline = -g_ca * ( m_infinity_v.* (v-v_ca) ) - g_l * (v-v_l) + 30;
    
    v_null_cline = numerator_v_null_cline./denominator_v_null_cline;

    % MLE 2nd set of variables, I ext = 30 quiver plot and  Null clines 
    figure(7)
        hold on
        plot(v, 100*w_null_cline)
        plot(v, 100*v_null_cline)

        % quiver plots
        [v_quiver,w_quiver] = meshgrid(linspace(-80,80,30), linspace(0,1,30));
        m_infinity_v_quiver = 0.5 * ( 1 + tanh((v_quiver-v1)/(v2)) ); 

        tau_w = 1./cosh((v_quiver-v3)/(2*v4));
        dv_dt = (1/c)*((-g_ca * ( m_infinity_v_quiver.*(v_quiver-v_ca) )) + (-g_k * ( w_quiver.*(v_quiver-v_k) )) + (-g_l * (v_quiver - v_l)) + 30);
        dw_dt = phi * (0.5 * ( 1 + tanh((v_quiver-v3)/(v4)) ) - w_quiver)./tau_w;
        quiver(v_quiver,100*w_quiver, dv_dt, 100*dw_dt, 1, 'color',[0 0 0]); % arrow length scaled 2 times for visibility

        % to draw manifolds
        % find eigen vectors at saddle node
        syms v_var3 w_var3
        dv_dt3 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var3-v1)/(v2)) ))*(v_var3-v_ca) )) + (-g_k * ( w_var3*(v_var3-v_k) )) + (-g_l * (v_var3 - v_l)) + 30);
        dw_dt3 = phi * (0.5 * ( 1 + tanh((v_var3-v3)/(v4)) ) - w_var3)/(1/cosh((v_var3-v3)/(2*v4)));
        
        df1_dv3 = diff(dv_dt3, v_var3);
        df1_dw3 = diff(dv_dt3, w_var3);
        df2_dv3 = diff(dw_dt3, v_var3);
        df2_dw3 = diff(dw_dt3, w_var3);
    
        v_eq3 = -19.5632
        w_eq3 = 0.0259
    
        % jacobian matrix and their eigen values
        jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
        [eigen_vec, eigen_val] = eig(jacobian3);
        disp("eigen vecs")
        disp(double(eigen_vec))
        disp(double(eigen_val))

         % [t r] = mle_solution_i_ext_set2(30, -19.56, 0.02)  
        % plot(r(:,1), 100*r(:,2))      

        hold off
    grid

    % finding the new equilibrium points
    new_eqs = []
    F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)) + 30);  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
    starting_pt = [-41; 0.02];
    options = optimoptions('fsolve','Display','iter');
    [x,fval] = fsolve(F,starting_pt,options);
    disp(x);
    new_eqs = [new_eqs, x];

    F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)) + 30);  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
    starting_pt = [-20; 2];
    options = optimoptions('fsolve','Display','iter');
    [x,fval] = fsolve(F,starting_pt,options);
    disp(x);
    new_eqs = [new_eqs, x];


    F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)) + 30);  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
    starting_pt = [4; 28];
    options = optimoptions('fsolve','Display','iter');
    [x,fval] = fsolve(F,starting_pt,options);
    disp(x);
    new_eqs = [new_eqs, x];
    
    
    
    disp("new equilibrium pts, new set of MLE, iext = 30")
    disp(new_eqs)
  % characterising the new equilbrium points
  for i=1:3
    disp("for point")
    % finding jacobian values
    syms v_var3 w_var3
    dv_dt3 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var3-v1)/(v2)) ))*(v_var3-v_ca) )) + (-g_k * ( w_var3*(v_var3-v_k) )) + (-g_l * (v_var3 - v_l)) + 30);
    dw_dt3 = phi * (0.5 * ( 1 + tanh((v_var3-v3)/(v4)) ) - w_var3)/(1/cosh((v_var3-v3)/(2*v4)));
    
    df1_dv3 = diff(dv_dt3, v_var3);
    df1_dw3 = diff(dv_dt3, w_var3);
    df2_dv3 = diff(dw_dt3, v_var3);
    df2_dw3 = diff(dw_dt3, w_var3);

    v_eq3 = new_eqs(1,i)
    w_eq3 = new_eqs(2,i)

    % jacobian matrix and their eigen values
    jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
    eigen_values3 = double(eig(jacobian3)) 
    disp("--------------------------")
end

    % frequency of action potential vs current
    figure(12)
        rates_of_ap = [];
        i_ext = [];
        for i=30:45
            i_ext = [i_ext, i];
            [t r] = mle_solution_i_ext_set2(i, 100, 0.10);
            frequency_of_ap = 1/calculate_ap_time(r,t);
            rates_of_ap = [rates_of_ap, frequency_of_ap];
        end

        plot(i_ext, rates_of_ap);
    grid


end



function ap_time = calculate_ap_time(r,t)
        rounded_off_voltages = round(r(:,1), 0);
        freq_table_all = tabulate(round(r(:,1), 0));
        
        [s, s1] = size(freq_table_all);
        % just remove the negative, we are calculating from positive
        freq_table = zeros(s,s1);
        for i=1:s
            if freq_table_all(i,1) > 0
                freq_table(i,1) = freq_table_all(i,1) ;
                freq_table(i,2) = freq_table_all(i,2) ;
                freq_table(i,3) = freq_table_all(i,3) ;
            end
        end
        
        % find the max_frequency maximum value
        maximum_frequency = freq_table(1,2);
        voltage_having_max_freq = freq_table(1,1);
        for i=1:s
            if freq_table(i,2) > maximum_frequency
                maximum_frequency = freq_table(i,2);
                voltage_having_max_freq = freq_table(i, 1);
            end
        end
        
        % find all the locations of max-positive_frequency
        locations_of_max_positive_frequency = [];
        [r_rows ,r_cols]= size(rounded_off_voltages);
        for i=1:r_rows
            if rounded_off_voltages(i, 1) == voltage_having_max_freq;
                locations_of_max_positive_frequency = [locations_of_max_positive_frequency, i];
            end
        end 

        % a peak in middle
        [l_row, l_col] = size(locations_of_max_positive_frequency);
        random_point_in_middle = floor(l_col/2);
        peak_in_middle_location = locations_of_max_positive_frequency(1,random_point_in_middle);
        
        
        t0 = t(peak_in_middle_location,1);
        negative_reached = 0;
        for i=peak_in_middle_location+1:r_rows
            if rounded_off_voltages(i,1) > rounded_off_voltages(i-1,1) & negative_reached
                break
            end

            if rounded_off_voltages(i,1) < 0
                negative_reached = 1;
            end
        end

        t1 = t(i,1);
       
        positive_reached = 0;
        for j=i+1:r_rows
            if rounded_off_voltages(j,1) < rounded_off_voltages(j-1) & positive_reached
                break
            end

            if rounded_off_voltages(j,1) > 0
                positive_reached = 1;
            end

        end

        t2 = t(j,1);
        
        ap_time = t2 - t0;

end


function result = mle_diff_eqn(t,r)

    % defining first set of MLE variables
    g_ca = 4.4;
    g_k = 8;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.02;
    v1 = -1.2;
    v2 = 18;
    v3 = 2;
    v4 = 30;
    v5 = 2;
    v6 = 30;
    c = 20;

    result = zeros(2,1);
    result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)));
    result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
end

function result = mle_diff_eqn_with_i_ext_steady(t,r)

    % defining first set of MLE variables
    g_ca = 4.4;
    g_k = 8;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.02;
    v1 = -1.2;
    v2 = 18;
    v3 = 2;
    v4 = 30;
    v5 = 2;
    v6 = 30;
    c = 20;

    result = zeros(2,1);
    result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)) + 86);
    result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
end

function result = mle_diff_eqn_with_i_ext_steady_backward_time(t,r)

    % defining first set of MLE variables
    g_ca = 4.4;
    g_k = 8;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.02;
    v1 = -1.2;
    v2 = 18;
    v3 = 2;
    v4 = 30;
    v5 = 2;
    v6 = 30;
    c = 20;

    result = zeros(2,1);
    result(1) = (-1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)) + 86);
    result(2) = -phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
end

function [t_vec,r_vec] = mle_solution_i_ext_set1(i_ext)
    
    
    function result = mle_diff_eqn_with_i_ext_steady2(t,r)

        % defining first set of MLE variables
        g_ca = 4.4;
        g_k = 8;
        g_l = 2;
        v_ca = 120;
        v_k = -84;
        v_l = -60;
        phi = 0.02;
        v1 = -1.2;
        v2 = 18;
        v3 = 2;
        v4 = 30;
        v5 = 2;
        v6 = 30;
        c = 20;
    
        result = zeros(2,1);
        result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)) + i_ext);
        result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
    end 
    
    [t_vec r_vec] = ode15s(@mle_diff_eqn_with_i_ext_steady2, [0 10000], [-80 0.10]);

end


function result = mle_diff_eqn2(t,r)

    % defining first set of MLE variables, a different phi
    g_ca = 4.4;
    g_k = 8;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.04;
    v1 = -1.2;
    v2 = 18;
    v3 = 2;
    v4 = 30;
    v5 = 2;
    v6 = 30;
    c = 20;

    result = zeros(2,1);
    result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)));
    result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
end

function result = mle_diff_eqn3(t,r)

    % defining first set of MLE variables, a different phi
    g_ca = 4.4;
    g_k = 8;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.01;
    v1 = -1.2;
    v2 = 18;
    v3 = 2;
    v4 = 30;
    v5 = 2;
    v6 = 30;
    c = 20;

    result = zeros(2,1);
    result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)));
    result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
end

function result = mle_diff_eqn_with_i_ext_steady_second_set(t,r)

    % defining second set of MLE variables
    g_ca = 4;
    g_k = 8.0;
    g_l = 2;
    v_ca = 120;
    v_k = -84;
    v_l = -60;
    phi = 0.0667;
    v1 = -1.2;
    v2 = 18;
    v3 = 12;
    v4 = 17.4;
    v5 = 12;
    v6 = 17.4;
    c = 20;

    result = zeros(2,1);
    result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)) + 30);
    result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
end


function [t_vec,r_vec] = mle_solution_i_ext_set2(i_ext, v_0, w_0)
    
    
    function result = mle_diff_eqn_with_i_ext_steady_second_set(t,r)

        % defining second set of MLE variables
        g_ca = 4;
        g_k = 8.0;
        g_l = 2;
        v_ca = 120;
        v_k = -84;
        v_l = -60;
        phi = 0.0667;
        v1 = -1.2;
        v2 = 18;
        v3 = 12;
        v4 = 17.4;
        v5 = 12;
        v6 = 17.4;
        c = 20;
    
        result = zeros(2,1);
        result(1) = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)) + i_ext);
        result(2) = phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
    end 
    
    [t_vec r_vec] = ode15s(@mle_diff_eqn_with_i_ext_steady_second_set, [0 10000], [v_0 w_0]);

end

function [t_vec,r_vec] = mle_solution_i_ext_set2_backward_time(i_ext, v_0, w_0)
    
    
    function result = mle_diff_eqn_with_i_ext_steady_second_set(t,r)

        % defining second set of MLE variables
        g_ca = 4;
        g_k = 8.0;
        g_l = 2;
        v_ca = 120;
        v_k = -84;
        v_l = -60;
        phi = 0.0667;
        v1 = -1.2;
        v2 = 18;
        v3 = 12;
        v4 = 17.4;
        v5 = 12;
        v6 = 17.4;
        c = 20;
    
        result = zeros(2,1);
        result(1) = (-1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((r(1)-v1)/(v2)) ))*(r(1)-v_ca) )) + (-g_k * ( r(2)*(r(1)-v_k) )) + (-g_l * (r(1) - v_l)) + i_ext);
        result(2) = -phi * (0.5 * ( 1 + tanh((r(1)-v3)/(v4)) ) - r(2))/(1/cosh((r(1)-v3)/(2*v4)));
    end 
    
    [t_vec r_vec] = ode15s(@mle_diff_eqn_with_i_ext_steady_second_set, [0 10000], [v_0 w_0]);

end