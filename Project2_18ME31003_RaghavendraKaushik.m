function project
    v = linspace(-80, 80);

    % defining variables
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

      
        % finding point of intersection
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
        
        [v_quiver,w_quiver] = meshgrid(linspace(-80,80,30), linspace(0,1,30));
        m_infinity_v_quiver = 0.5 * ( 1 + tanh((v_quiver-v1)/(v2)) ); 

        tau_w = 1./cosh((v_quiver-v3)/(2*v4));
        dv_dt = (1/c)*((-g_ca * ( m_infinity_v_quiver.*(v_quiver-v_ca) )) + (-g_k * ( w_quiver.*(v_quiver-v_k) )) + (-g_l * (v_quiver - v_l)));
        dw_dt = phi * (0.5 * ( 1 + tanh((v_quiver-v3)/(v4)) ) - w_quiver)./tau_w;
        quiver(v_quiver,100*w_quiver, dv_dt, 100*dw_dt, 7, 'color',[0 0 0]); % arrow length scaled 2 times for visibility
        
        % defininf dv/dt and dw/dt again because this time we want expressions not vectors for the jacobian
        % dv/dt and dw_dt are functions of (v_var, w_var)
        syms v_var w_var
        dv_dt2 = (1/c)*((-g_ca * ( (0.5 * ( 1 + tanh((v_var-v1)/(v2)) ))*(v_var-v_ca) )) + (-g_k * ( w_var*(v_var-v_k) )) + (-g_l * (v_var - v_l)));
        dw_dt2 = phi * (0.5 * ( 1 + tanh((v_var-v3)/(v4)) ) - w_var)/(1/cosh((v_var-v3)/(2*v4)));
        
        df1_dv = diff(dv_dt2, v_var);
        df1_dw = diff(dv_dt2, w_var);
        df2_dv = diff(dw_dt2, v_var);
        df2_dw = diff(dw_dt2, w_var);

        
        % jacobian matrix
        jacobian = [subs(df1_dv,{v_var,w_var},{v_eq, w_eq}) subs(df1_dw,{v_var,w_var},{v_eq, w_eq}); subs(df2_dv,{v_var,w_var},{v_eq, w_eq}) subs(df2_dw,{v_var,w_var},{v_eq, w_eq})  ]
        eigen_values = double(eig(jacobian))


    hold off

end

