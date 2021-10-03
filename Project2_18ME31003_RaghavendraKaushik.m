

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
        
        % running the equation from Equilibrium point
        [t,r] = ode15s(@mle_diff_eqn,[0 50000],[x(1) x(2)])
        plot(r(:,1),r(:,2))

    hold off

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