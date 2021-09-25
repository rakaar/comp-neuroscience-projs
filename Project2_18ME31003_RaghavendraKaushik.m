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

        plot(v, w_null_cline)
        plot(v, v_null_cline)

      
        % finding point of intersection
        for cnt=1:numel(v)
            s1 = v_null_cline(cnt);
            s2 = w_null_cline(cnt);
            zz = abs(s1-s2);
            if zz <= 0.01
               fprintf("v_eq,w_eq %.3f, %.3f \n", v(cnt), s1)
            end
        end

    hold off
end