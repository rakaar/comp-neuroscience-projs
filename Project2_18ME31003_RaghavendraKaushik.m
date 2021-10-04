

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
    
    
    figure(1)
    % plotting null-clines,finding equilibrium points and their nature, quiver plots
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
        jacobian = [subs(df1_dv,{v_var,w_var},{v_eq, w_eq}) subs(df1_dw,{v_var,w_var},{v_eq, w_eq}); subs(df2_dv,{v_var,w_var},{v_eq, w_eq}) subs(df2_dw,{v_var,w_var},{v_eq, w_eq})  ];
        eigen_values = double(eig(jacobian)) % we see that eigen values are negative, implying that equilibrium point is a stable point
        
        % running the equation from Equilibrium point, it stays there
        [t,r] = ode15s(@mle_diff_eqn,[0 100],[x(1) x(2)])
        plot(r(:,1),100*r(:,2))

    hold off
    grid

    % plotting for different values of phi
    figure(2)
        hold on
        % [t1,r1] = ode15s(@mle_diff_eqn,[0 10000],[x(1)+10 x(2)]) % sub threshold, phi = 0.02
        % plot4 = plot(r1(:,1),r1(:,2))

        

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

        
        % Calculations verifying the plots - finding new equilibrium point, its jacobian's eigen values to know its kind
        % estimating the point of intersection with iext = 86
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
        hold off
    grid


    % other all i currents comparison
    iext = [80, 86, 90]
    for i=1:3
        disp("for i ext is ")
        disp(iext(i))
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

    figure(7)
        plot(v, 100*w_null_cline)
        plot(v, 100*v_null_cline)
    grid



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