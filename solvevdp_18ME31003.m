% A. Raghavendra Kaushik - 18ME31003
function project
    % tspan = [0 1000]; % for phase plane and time comparison
    tspan = [0 50]; % for plotting indvidual solutions, so that figure isn't clumsy
    y0 = [1 0];

    disp('mu = 0.1')
    tic, [t,v] = ode45(@vdp_pt_1,tspan,y0);
    toc


    disp('mu = 1')
    tic, [t1,v1] = ode45(@vdp_1,tspan,y0);
    toc

    disp('mu = 100 - ode45')
    tic, [t2,v2] = ode45(@vdp_100,tspan,y0);
    toc

    disp('mu = 100 - ode15s')
    tic, [t2s,v2s] = ode15s(@vdp_100,tspan,y0);
    toc

    % plotting y(t)
    figure(1)
        subplot(1,3,1)
        plot(t,v(:,1),t,v(:,2))
        title('mu=0.1')

        subplot(1,3,2)
        plot(t1,v1(:,1),t1,v1(:,2))
        title('mu=1')

        subplot(1,3,3)
        plot(t2,v2(:,1),t2,v2(:,2))
        title('mu=100')
    grid


    % phase plane [y,z]
    figure(2)
        subplot(1,3,1)
        plot(v(:,1),v(:,2))
        title('mu=0.1')

        subplot(1,3,2)
        plot(v1(:,1),v1(:,2))
        title('mu=1')

        subplot(1,3,3)
        plot(v2s(:,1),v2s(:,2))
        title('mu=100')
    grid

end

% defining the vdp eqns with different mu values 
% y' = f(t), z'=f(t) is the system. Returns [y z]'
function result = vdp_pt_1(t,r)
    mu = 0.1;
    result = zeros(2,1);
    result(1) = mu*r(2) ;% y
    result(2) = mu*r(2)*(1-r(1)^2) - r(1)/mu ; % z
end

function result = vdp_1(t,r)
    mu = 1;
    result = zeros(2,1);
    result(1) = mu*r(2) ;% y
    result(2) = mu*r(2)*(1-r(1)^2) - r(1)/mu ; % z
end

function result = vdp_100(t,r)
    mu = 100;
    result = zeros(2,1);
    result(1) = mu*r(2) ;% y
    result(2) = mu*r(2)*(1-r(1)^2) - r(1)/mu ; % z
end

