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
            for step_size = 45.5:0.01:46.5
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
        quiver(v_quiver,100*w_quiver, dv_dt, 100*dw_dt, 3, 'color',[0 0 0]); % arrow length scaled 3 times for visibility

        % to draw manifolds
        % find eigen vectors at saddle node - unstable manifolds
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>> forward in time >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
      

        length_of_eigen_vec1 = (double(eigen_vec(1,1))^2 + double(eigen_vec(2,1))^2)^0.5;
        x_pt_on_eigen_vec_1_nearer_to_saddle_node = v_eq3 + 0.1*(double(eigen_vec(1,1))/length_of_eigen_vec1);
        y_pt_on_eigen_vec_1_nearer_to_saddle_node = w_eq3 + 0.1*(double(eigen_vec(2,1))/length_of_eigen_vec1);
        
        
        [t3 r3] = mle_solution_i_ext_set2(30, x_pt_on_eigen_vec_1_nearer_to_saddle_node, y_pt_on_eigen_vec_1_nearer_to_saddle_node);
        unstableManifold1 = plot(r3(:,1), 100*r3(:,2));
        
        
        length_of_eigen_vec2 = (double(eigen_vec(1,2))^2 + double(eigen_vec(2,2))^2)^0.5;
        x_pt_on_eigen_vec_2_nearer_to_saddle_node = v_eq3 + 0.2*(double(eigen_vec(1,2))/length_of_eigen_vec2);
        y_pt_on_eigen_vec_2_nearer_to_saddle_node = w_eq3 + 0.2*(double(eigen_vec(2,2))/length_of_eigen_vec2);
        [t4 r4] = mle_solution_i_ext_set2(30, x_pt_on_eigen_vec_2_nearer_to_saddle_node, y_pt_on_eigen_vec_2_nearer_to_saddle_node);
        unstableManifold2 = plot(r4(:,1), 100*r4(:,2));

        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< backward in time <<<<<<<<<<<<<<<<<<<<<<<
        disp("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<backward in time <<<<<<<<<<<<<<<<")
        % a different point to start
        x_pt_on_eigen_vec_1_nearer_to_saddle_node = v_eq3 - 0.002*(double(eigen_vec(1,1))/length_of_eigen_vec1);
        y_pt_on_eigen_vec_1_nearer_to_saddle_node = w_eq3 - 0.002*(double(eigen_vec(2,1))/length_of_eigen_vec1);
        [t1 r1] = mle_solution_i_ext_set2_backward_time(30,x_pt_on_eigen_vec_1_nearer_to_saddle_node, y_pt_on_eigen_vec_1_nearer_to_saddle_node);
        stable_manifold1 = plot(r1(:,1), 100*r1(:,2));

        [t2 r2] = mle_solution_i_ext_set2_backward_time(30,  x_pt_on_eigen_vec_2_nearer_to_saddle_node, y_pt_on_eigen_vec_2_nearer_to_saddle_node);
        stable_manifold2 = plot(r2(:,1), 100*r2(:,2));

        legend([stable_manifold1, stable_manifold2, unstableManifold1, unstableManifold2], ["stable manifold 1", "stable manifold 2", "unstable manifold 1", "unstable manifold 2"]);
        

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

    % finding equilibrium points for range 30 to 50

    disp("analysing equilbirum points from 30 to 50");
    for i=30:50
        fprintf("i ext is %f \n",i);
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

            syms v_var3 w_var3
            dv_dt3 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var3-v1)/(v2)) ))*(v_var3-v_ca) )) + (-g_k * ( w_var3*(v_var3-v_k) )) + (-g_l * (v_var3 - v_l)) + i);
            dw_dt3 = phi * (0.5 * ( 1 + tanh((v_var3-v3)/(v4)) ) - w_var3)/(1/cosh((v_var3-v3)/(2*v4)));
            
            df1_dv3 = diff(dv_dt3, v_var3);
            df1_dw3 = diff(dv_dt3, w_var3);
            df2_dv3 = diff(dw_dt3, v_var3);
            df2_dw3 = diff(dw_dt3, w_var3);
        
            F = @(x) [(1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((x(1)-v1)/(v2)) ))*(x(1)-v_ca) )) + (-g_k * ( x(2)*(x(1)-v_k) )) + (-g_l * (x(1) - v_l)) + i);  phi * (0.5 * ( 1 + tanh((x(1)-v3)/(v4)) ) - x(2))/(1/cosh((x(1)-v3)/(2*v4)))];
            options = optimoptions('fsolve','Display','off');
            if i < 40 
                starting_pt = [-41; 0.02];
                [x,fval] = fsolve(F,starting_pt,options);
                v_eq3 = x(1);
                w_eq3 = x(2);
                jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
                eigen_values3 = double(eig(jacobian3));
                disp(x);
                disp(eigen_values3)
                    
                starting_pt = [-20; 2];
                [x,fval] = fsolve(F,starting_pt,options);
                v_eq3 = x(1);
                w_eq3 = x(2);
                disp(x);
                jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
                eigen_values3 = double(eig(jacobian3));
                disp(eigen_values3)
            end
            starting_pt = [4; 28];
            [x,fval] = fsolve(F,starting_pt,options);
            v_eq3 = x(1);
            w_eq3 = x(2);
            disp(x);
            jacobian3 = [subs(df1_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df1_dw3,{v_var3,w_var3},{v_eq3, w_eq3}); subs(df2_dv3,{v_var3,w_var3},{v_eq3, w_eq3}) subs(df2_dw3,{v_var3,w_var3},{v_eq3, w_eq3})  ];
            eigen_values3 = double(eig(jacobian3));
            disp(eigen_values3)
                
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hogkin Huxley %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find out value of E leak such that resting potential is zero
    % mathematically, dv/dt = 0, at v = -60 and E leak = ? and some values of m, n, h
    % m, n, h are found by dm/dt = 0, dn/dt = 0; dh/dt = 0
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hogkin Huxley %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5");
    alpha_m = (-0.1 * (-60 + 35))/(exp(-(-60 + 35)/10) - 1);
    beta_m = 4;

    alpha_n = (0.1)/(exp(1) - 1);
    beta_n = 0.125; 

    alpha_h = 0.07 ;
    beta_h = 1/(1 + exp(3));

    m_eq = alpha_m/(alpha_m + beta_m);
    h_eq = alpha_h/(alpha_h + beta_h);
    n_eq = alpha_n/(alpha_n + beta_n);

    fprintf("m h n are %f %f %f \n", m_eq, h_eq, n_eq);

     % vars
     g_k_bar = 36;
     e_k = -72;
 
     g_na_bar = 120;
     e_na = 55;
 
     g_l = 0.3;
     
 
     c = 1;
     
     v_eq = -60;

     % calculating e_leak for i ext = 0
     e_leak = v_eq + ((1/g_l) * ((g_k_bar * n_eq^4 * (v_eq - e_k)) + (g_na_bar * m_eq^3 * h_eq * (v_eq - e_na))));
    fprintf("e_leak is %f \n", e_leak);
    % e_leak = -49.401079
    
    figure(13)
        disp("plotting hh")
        [t_h,r_h] = ode15s(@hh,[0 1000],[10  0.052932 0.596121 0.317677]);
        plot(t_h, r_h(:,1));
    grid


    % finding out stability of hh model, i ext = 0
    % 1. find the point, 2. find eigen values
    % x = [v m h n]
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;
    
    % solving equilbirum point using vpasolve 
    syms v m h n
    X = vpasolve([
        (g_k_bar * (n^4) * (v - e_k))  + (g_na_bar * (m^3) * h * (v - e_na)) + (g_l * (v - e_l)) == 0,
        ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0,
        (0.07 * exp(-(v+60)/20))*(1-h) - (1/(exp(-(v+30)/10) + 1))*(h) == 0,
        ((-0.01 * (v+50))/(exp(-(v+50)/10) - 1))*(1-n) - (0.125 * (exp(-(v+60)/80)))*(n) == 0
    ], [v,m,h,n]);
    disp("equilibrium points for hh ");
    fprintf("v m h n %f %f %f %f \n", X.v, X.m, X.h, X.n);

    % finding eigen values of jacobian
    syms v1 m1 h1 n1
    dv_dt = (1/c)* (-(g_k_bar * (n1^4) * (v1 - e_k))  - (g_na_bar * (m1^3) * h1 * (v1 - e_na)) - (g_l * (v1 - e_l)));
    dm_dt = ((-0.1 * (v1+35))/(exp(-(v1+35)/10) -1))*(1-m1) - (4 * exp(-(v1+60)/18))*(m1);
    dh_dt = (0.07 * exp(-(v1+60)/20))*(1-h1) - (1/(exp(-(v1+30)/10) + 1))*(h1);
    dn_dt = ((-0.01 * (v1+50))/(exp(-(v1+50)/10) - 1))*(1-n1) - (0.125 * (exp(-(v1+60)/80)))*(n1);

    jacobian = zeros(4,4);
    jacobian(1,1) = subs(diff(dv_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(1,2) = subs(diff(dv_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(1,3) = subs(diff(dv_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(1,4) = subs(diff(dv_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    jacobian(2,1) = subs(diff(dm_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(2,2) = subs(diff(dm_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(2,3) = subs(diff(dm_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(2,4) = subs(diff(dm_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    jacobian(3,1) = subs(diff(dh_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(3,2) = subs(diff(dh_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(3,3) = subs(diff(dh_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(3,4) = subs(diff(dh_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    jacobian(4,1) = subs(diff(dn_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(4,2) = subs(diff(dn_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(4,3) = subs(diff(dn_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(4,4) = subs(diff(dn_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    eigen_values = double(eig(jacobian));
    disp("eigen values of hh model 4d");
    disp(eigen_values);

    % % calculating threshold for hh
    % Non zero AP not a good idea, just see v rise vs time for intial few seconds
    figure(141)
        v_max = [];
        v_inital_values = [];
        for step_size=1:20
            [t r] = ode15s(@hh, [0 1000] ,[-60+step_size  0.052932 0.596121 0.317677]);
            v_max = [v_max, max(r(:,1))];
            v_inital_values = [v_inital_values, -60+step_size];
        end
        plot(v_inital_values, v_max)
    grid

   
    % stability of equilbirum points
    for i=8:12
        syms v m h n
    X = vpasolve([
        -i + (g_k_bar * (n^4) * (v - e_k))  + (g_na_bar * (m^3) * h * (v - e_na)) + (g_l * (v - e_l)) == 0,
        ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0,
        (0.07 * exp(-(v+60)/20))*(1-h) - (1/(exp(-(v+30)/10) + 1))*(h) == 0,
        ((-0.01 * (v+50))/(exp(-(v+50)/10) - 1))*(1-n) - (0.125 * (exp(-(v+60)/80)))*(n) == 0
    ], [v,m,h,n]);
    fprintf("I ext is %f\n",i);
    fprintf("v m h n %f %f %f %f \n", X.v, X.m, X.h, X.n);

    % finding eigen values of jacobian
    syms v1 m1 h1 n1
    dv_dt = (1/c)* (i-(g_k_bar * (n1^4) * (v1 - e_k))  - (g_na_bar * (m1^3) * h1 * (v1 - e_na)) - (g_l * (v1 - e_l)));
    dm_dt = ((-0.1 * (v1+35))/(exp(-(v1+35)/10) -1))*(1-m1) - (4 * exp(-(v1+60)/18))*(m1);
    dh_dt = (0.07 * exp(-(v1+60)/20))*(1-h1) - (1/(exp(-(v1+30)/10) + 1))*(h1);
    dn_dt = ((-0.01 * (v1+50))/(exp(-(v1+50)/10) - 1))*(1-n1) - (0.125 * (exp(-(v1+60)/80)))*(n1);

    jacobian = zeros(4,4);
    jacobian(1,1) = subs(diff(dv_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(1,2) = subs(diff(dv_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(1,3) = subs(diff(dv_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(1,4) = subs(diff(dv_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    jacobian(2,1) = subs(diff(dm_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(2,2) = subs(diff(dm_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(2,3) = subs(diff(dm_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(2,4) = subs(diff(dm_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    jacobian(3,1) = subs(diff(dh_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(3,2) = subs(diff(dh_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(3,3) = subs(diff(dh_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(3,4) = subs(diff(dh_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    jacobian(4,1) = subs(diff(dn_dt, v1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(4,2) = subs(diff(dn_dt, m1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(4,3) = subs(diff(dn_dt, h1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});
    jacobian(4,4) = subs(diff(dn_dt, n1), {v1, m1, h1, n1}, {X.v, X.m, X.h, X.n});

    eigen_values = double(eig(jacobian));
    disp(eigen_values);
    end


    figure(16)
        hold on
            disp("Myotonic hh")
            
            subplot(2,2,1)
            [t r] = myotonoic_hh(0);
            plot(t, r(:,1));
            title('f_in = 0')

            subplot(2,2,2)
            [t1 r1] = myotonoic_hh(0.1);
            plot(t1, r1(:,1));
            title('f_in = 0.1');

            subplot(2,2,3)
            [t2 r2] = myotonoic_hh(0.17);
            plot(t2, r2(:,1));
            title('f_in = 0.17')

            subplot(2,2,4)
            [t3 r3] = myotonoic_hh(0.2);
            plot(t3, r3(:,1));
            title('f_in = 0.20')
        hold off
    grid

  
    % v-n reduced model
    figure(17)
        v_initals = [];
        v_maxs = [];
        for i=1:10
            [t r] = hh_2d_iext(i);
            v_initals = [v_initals, i];
            v_maxs = [v_maxs, max(r(:,1))];
        end
    plot(v_initals, v_maxs);
    grid

    figure(171)
    hold on
        for i=1:5
            [t r] = hh_2d_iext(i);
            plot(r(:,1), 100*r(:,2));
        end
    hold off
    grid



    
    % phase plane analysis n-v with myotonia
    
    figure(18)
        v18 = linspace(-70, 70);
            
        if v18 == -50 
            alpha_n = 0.1;
        else
            alpha_n = (-0.01 * (v18 + 50))./(exp(-(v18 + 50)/10) - 1);
        end
        beta_n = 0.125 * exp(-(v18 + 60)/80);

        if v18 == -35
            alpha_m = 1;
        else
            alpha_m = (-0.1 * (v18 + 35))./(exp(-(v18 + 35)/10) - 1);
        end
        beta_m = 4 * exp(-(v18 + 60)/18);
    
        m_inf = alpha_m./(alpha_m + beta_m);
        h = 0.596121;
    
        f_ni = 0.20;
        iext = 0;

        g_k_bar = 36;   e_k = -72;    g_na_bar = 120;    e_na = 55;    g_l = 0.3;   e_l = -49.401079;
       
        % 18_nullcline = (1/c) * (iext - (g_k_bar * (r(2)^4) * (v - e_k))   - (g_na_bar * (1-f_ni)* (m^3) * h * (r(1) - e_na)) - (g_na_bar * (f_ni)* (m^3) * (r(1) - e_na)) - (g_l * (r(1) - e_l)) ) ;
        v_null_cline18 = ((-(g_na_bar *(1-f_ni) * h * (m_inf.^3) .* (v18 - e_na)) - (g_na_bar * f_ni * (m_inf.^3) .* (v18 - e_na)) - (g_l * (v18 - e_l))) ./ (g_k_bar * (v18 - e_k))).^ (1/4) ;
        n_nullcline18 = alpha_n ./ (alpha_n + beta_n);
        hold on
            plot(v18, 100*v_null_cline18);
            plot(v18, 100*n_nullcline18);
        hold off

        initial_guess = zeros(3,2);

        initial_guess(1,1) = -60.1;
        initial_guess(1,2) = 0.31;


        initial_guess(2,1) = -47.3;
        initial_guess(2,2) = 0.51;


        initial_guess(3,1) = 19.09;
        initial_guess(3,2) = 0.93;


        for f_ni=0.02:+0.02:0.4
            fprintf("f_ni is %f  ",f_ni);
            for i=1:3
                syms v n
                X = vpasolve([
                     (g_k_bar * (n^4) * (v - e_k))  + (g_na_bar * (1-f_ni) * ((1/(1+((4 * exp(-(v+60)/18))/((-0.1 * (v + 35))/(exp(-(v + 35)/10) - 1)))))^3) * h * (v - e_na)) + (g_na_bar * (f_ni) * ((1/(1+((4 * exp(-(v+60)/18))/((-0.1 * (v + 35))/(exp(-(v + 35)/10) - 1)))))^3)  * (v - e_na)) + (g_l * (v - e_l)) == 0,
                    ((-0.01 * (v+50))/(exp(-(v+50)/10) - 1))*(1-n) - (0.125 * (exp(-(v+60)/80)))*(n) == 0
                ], [v, n], [initial_guess(i,1);initial_guess(i,2)]);
    
                % fprintf("v n  %f %f \n", X.v, X.n);
                
                syms  v1 n1
                dv_dt = (1/c)* (-(g_k_bar * (n1^4) * (v1 - e_k))  - (g_na_bar * (1-f_ni) * ((1/(1+((4 * exp(-(v1+60)/18))/((-0.1 * (v1 + 35))/(exp(-(v1 + 35)/10) - 1)))))^3) * h * (v1 - e_na)) - (g_na_bar * f_ni * ((1/(1+((4 * exp(-(v1+60)/18))/((-0.1 * (v1 + 35))/(exp(-(v1 + 35)/10) - 1)))))^3) * (v1 - e_na)) - (g_l * (v1 - e_l)));
                dn_dt = ((-0.01 * (v1+50))/(exp(-(v1+50)/10) - 1))*(1-n1) - (0.125 * (exp(-(v1+60)/80)))*(n1);
                
                jacobian = zeros(2,2);
                jacobian(1,1) = subs(diff(dv_dt, v1), {v1,n1}, {X.v,  X.n});
                jacobian(1,2) = subs(diff(dv_dt, n1), {v1, n1}, {X.v,  X.n});
                jacobian(2,1) = subs(diff(dn_dt, v1), {v1, n1}, {X.v,  X.n});
                jacobian(2,2) = subs(diff(dn_dt, n1), {v1, n1}, {X.v,  X.n});
                
                eigen_values = double(eig(jacobian));
                % disp(eigen_values);
                fprintf("  %s ", get_stability(eigen_values));
            end
            fprintf("\n");
        end

    grid

  

  % anode break
    figure(19)
        % vars
        g_k_bar = 36;
        e_k = -72;

        g_na_bar = 120;
        e_na = 55;

        g_l = 0.3;
        e_l = -49.401079;
        c = 1;
        syms v m h n
        X = vpasolve([
            (g_k_bar * (n^4) * (v - e_k))  + (g_na_bar * (m^3) * h * (v - e_na)) + (g_l * (v - e_l)) == 0,
            ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0,
            (0.07 * exp(-(v+60)/20))*(1-h) - (1/(exp(-(v+30)/10) + 1))*(h) == 0,
            ((-0.01 * (v+50))/(exp(-(v+50)/10) - 1))*(1-n) - (0.125 * (exp(-(v+60)/80)))*(n) == 0
        ], [v,m,h,n]);

        v_rest = double(X.v);
        [t r] = ode15s(@hh_i_negative_for_anode_break, [0 50], [double(X.v) double(X.m) double(X.h) double(X.n)]);
        plot(t, r(:,1));
    grid



    syms v m h n
    X = vpasolve([
        -3 + (-g_k_bar * (n^4) * (v - e_k))  + (-g_na_bar * (m^3) * h * (v - e_na)) + (-g_l * (v - e_l)) == 0,
        ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0,
        (0.07 * exp(-(v+60)/20))*(1-h) - (1/(exp(-(v+30)/10) + 1))*(h) == 0,
        ((-0.01 * (v+50))/(exp(-(v+50)/10) - 1))*(1-n) - (0.125 * (exp(-(v+60)/80)))*(n) == 0
    ], [v,m,h,n]);

    v_h = double(X.v);
    

  
    % case 1
    v20 = linspace(-72,55);

    if v20 == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (v20 + 35))./(exp(-(v20 + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(v20+60)/18);
    m_nullcline = alpha_m./(alpha_m + beta_m);

    figure(203)
            n1 = get_n(v_rest);
            h1 = get_h(v_rest);

            v_nullclinef = @(V) 100*((((-g_k_bar * (n1^4) * (V - e_k)) + (-g_l * (V - e_l)))/(g_na_bar * h1 * (V - e_na)))^(1/3));
            hold on 
                fplot(@(V) v_nullclinef(V), [-72 55]);
                xlim([-72, 55]);
                ylim([0, 100]);

                plot(v20, 100*m_nullcline);
                xlim([-72, 55]);
                ylim([0, 100]);
            hold off


            initial_guess = zeros(3,2);

            initial_guess(1,1) = -60.1;
            initial_guess(1,2) = 0.05;

            initial_guess(2,1) = -57.1;
            initial_guess(2,2) = 0.07;

            initial_guess(3,1) = 53;
            initial_guess(3,2) = 0.98;

            for i=1:3
                syms v m
                X = vpasolve([
                    (-g_k_bar * (n1^4) * (v - e_k))  + (-g_na_bar * (m^3) * h1 * (v - e_na)) + (-g_l * (v - e_l)) == 0,
                    ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0
                ], [v, m], [initial_guess(i,1);initial_guess(i,2)]);

                
                
                syms  v1 m1
                dv_dt = (1/c)* (-(g_k_bar * (n1^4) * (v1 - e_k))  - (g_na_bar * (m1^3) * h1 * (v1 - e_na)) - (g_l * (v1 - e_l)));
                dm_dt = ((-0.1 * (v1+35))/(exp(-(v1+35)/10) -1))*(1-m1) - (4 * exp(-(v1+60)/18))*(m1);
                
                jacobian = zeros(2,2);
                jacobian(1,1) = subs(diff(dv_dt, v1), {v1,m1}, {X.v,  X.m});
                jacobian(1,2) = subs(diff(dv_dt, m1), {v1, m1}, {X.v,  X.m});
                jacobian(2,1) = subs(diff(dm_dt, v1), {v1, m1}, {X.v,  X.m});
                jacobian(2,2) = subs(diff(dm_dt, m1), {v1, m1}, {X.v,  X.m});
                
                eigen_values = double(eig(jacobian));
                
                fprintf("\n eigen values %f %f %s \n",eigen_values(1,1), eigen_values(2,1), get_stability(eigen_values));
            end
            
    grid
    
  


    figure(204)
    n2 = get_n(v_h);
    h2 = get_h(v_h);

    v_nullclinef2 = @(V) 100*((((-g_k_bar * (n2^4) * (V - e_k)) + (-g_l * (V - e_l)))/(g_na_bar * h2 * (V - e_na)))^(1/3));
    hold on 
        fplot(@(V) v_nullclinef2(V), [-72 55]);
        xlim([-72, 55]);
        ylim([0, 100]);
       
        plot(v20, 100*m_nullcline);
        xlim([-72, 55]);
        ylim([0, 100]);
    hold off
    grid

    syms v m
    X = vpasolve([
        (-g_k_bar * (n2^4) * (v - e_k))  + (-g_na_bar * (m^3) * h2 * (v - e_na)) + (-g_l * (v - e_l)) == 0,
        ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0
    ], [v, m], [53;0.98]);

    
    
    syms  v1 m1
    dv_dt = (1/c)* (-(g_k_bar * (n2^4) * (v1 - e_k))  - (g_na_bar * (m1^3) * h2 * (v1 - e_na)) - (g_l * (v1 - e_l)));
    dm_dt = ((-0.1 * (v1+35))/(exp(-(v1+35)/10) -1))*(1-m1) - (4 * exp(-(v1+60)/18))*(m1);
    
    jacobian = zeros(2,2);
    jacobian(1,1) = subs(diff(dv_dt, v1), {v1,m1}, {X.v,  X.m});
    jacobian(1,2) = subs(diff(dv_dt, m1), {v1, m1}, {X.v,  X.m});
    jacobian(2,1) = subs(diff(dm_dt, v1), {v1, m1}, {X.v,  X.m});
    jacobian(2,2) = subs(diff(dm_dt, m1), {v1, m1}, {X.v,  X.m});
    
    eigen_values = double(eig(jacobian));
    fprintf("\n eigen values %f %f %s \n",eigen_values(1,1), eigen_values(2,1), get_stability(eigen_values));

end

function n20 = get_n(v)
    if v == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (v + 50))/(exp(-(v + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(v + 60)/80);

    n20 = alpha_n/(alpha_n + beta_n);
end

function h20 = get_h(v)
    alpha_h = 0.07 * exp(-(v + 60)/20);
    beta_h = 1/(1 + exp(-(v+30)/10));
    
    h20 = alpha_h/(alpha_h + beta_h);
end

function stability_status = get_stability(eigen_values)
    v1 = eigen_values(1,1);
    v2 = eigen_values(2,1);

    if isreal(v1)
        if v1*v2 > 0
            if(v1 > 0)
                stability_status = "unstable";
            else
                stability_status = "stable";
            end
        else
            stability_status = "saddle";
        end

    else 
        if real(v1) > 0
            stability_status = "unstable spiral";
        else
            stability_status = "stable spiral";
        end
    end
end

function m_value = get_m(v_value)
    beta_m  = 4 * exp(-(v_value+60)/18);


    if v_value == -35
        m_value = 1/(1+beta_m);
    else
        alpha_m = (-0.1 * (v_value + 35))/(exp(-(v_value + 35)/10) - 1);

        m_value =  alpha_m/(alpha_m + beta_m);
    end
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
    
    
    function result = mle_diff_eqn_with_i_ext_steady_second_set_backward(t,r)

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
    
    [t_vec r_vec] = ode15s(@mle_diff_eqn_with_i_ext_steady_second_set_backward, [0 -1000], [v_0 w_0]);

end

%%%%%%%%%%%%%%%%% Hogkin Huxley Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = hh_2d(t,r)

    % vars
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;
    iext = 0;

   
    if r(1) == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(r(1) + 60)/80);

    if r(1) == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(r(1) + 60)/18);

    % values of m and h at equilibrium
    % m = 0.052932;
    h = 0.596121;

    
    result = zeros(2,1); % v,n
    result(1) = (1/c) * (iext - (g_k_bar * (r(2)^4) * (r(1) - e_k))   - (g_na_bar * ((alpha_m/(alpha_m + beta_m))^3) * h * (r(1) - e_na)) - (g_l * (r(1) - e_l)) ) ;
    result(2) = (alpha_n * (1 - r(2))) - (beta_n * r(2));
end


function [t_vec, r_vec] = myotonoic_hh_2d(f_ni)
    function result = hh_2d_defective(t,r)

        % vars
        g_k_bar = 36;
        e_k = -72;
    
        g_na_bar = 120;
        e_na = 55;
    
        g_l = 0.3;
        e_l = -49.401079;
    
        c = 1;
        iext = 10;
    
       
        if r(1) == -50 
            alpha_n = 0.1;
        else
            alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
        end
        beta_n = 0.125 * exp(-(r(1) + 60)/80);
    
        % values of m and h at equilibrium
        m = 0.052932;
        h = 0.596121;
    
        
        result = zeros(2,1); % v,n
        result(1) = (1/c) * (iext - (g_k_bar * (r(2)^4) * (r(1) - e_k))   - (g_na_bar * (1-f_ni)* (m^3) * h * (r(1) - e_na)) - (g_na_bar * (f_ni)* (m^3) * (r(1) - e_na)) - (g_l * (r(1) - e_l)) ) ;
        result(2) = (alpha_n * (1 - r(2))) - (beta_n * r(2));
    end 

    [t_vec r_vec] = ode15s(@hh_2d_defective, [0 300], [50  0.317677]);

end

function [t_vec,r_vec] = myotonoic_hh(f_ni)
    
    
    function result = hh_with_problem(t,r)

        % vars
        g_k_bar = 36;
        e_k = -72;

        g_na_bar = 120;
        e_na = 55;

        g_l = 0.3;
        e_l = -49.401079;

        c = 1;
        iext = 10;

        if r(1) == -35
            alpha_m = 1;
        else
            alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
        end
        beta_m = 4 * exp(-(r(1) + 60)/18);


        if r(1) == -50 
            alpha_n = 0.1;
        else
            alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
        end
        beta_n = 0.125 * exp(-(r(1) + 60)/80);

        alpha_h = 0.07 * exp(-(r(1) + 60)/20);
        beta_h = 1/(1 + exp(-(r(1)+30)/10));
        
        result = zeros(4,1); % v,m,h,n
        result(1) = (1/c) * ( iext - (g_k_bar * r(4)^4 * (r(1) - e_k)) - (g_na_bar * (1-f_ni) * r(2)^3 * r(3) * (r(1) - e_na)   +    g_na_bar * f_ni  * r(2)^3 *  (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
        result(2) = (alpha_m * (1 - r(2))) - (beta_m * r(2));   
        result(3) = (alpha_h * (1 - r(3))) - (beta_h * r(3));
        result(4) = (alpha_n * (1 - r(4))) - (beta_n * r(4));
    end 
    
    [t_vec r_vec] = ode15s(@hh_with_problem, [0 300], [50  0.052932 0.596121 0.317677]);

end

function [t_vec r_vec] = hh_2d_iext(v_step)
    function result = hh_2d_diffeqn(t,r)
        g_k_bar = 36;  e_k = -72; g_na_bar = 120; e_na = 55; g_l = 0.3;  e_l = -49.401079; c = 1;

        
        if r(1) == -35
            alpha_m = 1;
        else
            alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
        end
        beta_m = 4 * exp(-(r(1) + 60)/18);
        m_inf = alpha_m/(alpha_m + beta_m);


        if r(1) == -50 
            alpha_n = 0.1;
        else
            alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
        end
        beta_n = 0.125 * exp(-(r(1) + 60)/80);

        h = 0.596121;

        result = zeros(2,1);
        result(1) = (1/c) * ( 0 - (g_k_bar * r(2)^4 * (r(1) - e_k)) - (g_na_bar * (m_inf^3) * h * (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
        result(2) = (alpha_n * (1 - r(2))) - (beta_n * r(2));
    end

    [t_vec r_vec] = ode15s(@hh_2d_diffeqn, [0 1000], [-60+v_step 0.3]);
end

function result = hh(t,r)

    % vars
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;
    iext = 0;

   

    if r(1) == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(r(1) + 60)/18);


    if r(1) == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(r(1) + 60)/80);

    alpha_h = 0.07 * exp(-(r(1) + 60)/20);
    beta_h = 1/(1 + exp(-(r(1)+30)/10));
    
    result = zeros(4,1); % v,m,h,n
    result(1) = (1/c) * ( iext - (g_k_bar * r(4)^4 * (r(1) - e_k)) - (g_na_bar * r(2)^3 * r(3) * (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
    result(2) = (alpha_m * (1 - r(2))) - (beta_m * r(2));   
    result(3) = (alpha_h * (1 - r(3))) - (beta_h * r(3));
    result(4) = (alpha_n * (1 - r(4))) - (beta_n * r(4));
end



function result = hh_i_negative_for_anode_break(t,r)

    % vars
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;
   
   

    if r(1) == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(r(1) + 60)/18);


    if r(1) == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(r(1) + 60)/80);

    alpha_h = 0.07 * exp(-(r(1) + 60)/20);
    beta_h = 1/(1 + exp(-(r(1)+30)/10));
    
    if t < 20
        iext = -3;
    else
        iext = 0;
    end

    result = zeros(4,1); % v,m,h,n
    result(1) = (1/c) * (iext - (g_k_bar * r(4)^4 * (r(1) - e_k)) - (g_na_bar * r(2)^3 * r(3) * (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
    result(2) = (alpha_m * (1 - r(2))) - (beta_m * r(2));   
    result(3) = (alpha_h * (1 - r(3))) - (beta_h * r(3));
    result(4) = (alpha_n * (1 - r(4))) - (beta_n * r(4));
end

