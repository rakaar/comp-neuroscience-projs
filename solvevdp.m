tspan = [0 20];
y0 = [1 0];
[t,v] = ode45(@vdp,tspan,y0)

%aliter - using the inbuilt vdp1 in matlab
% [t,v] = ode45(@vdp1,tspan,y0)


plot(t,v(:,1),t,v(:,2))
% plot(v(:,1),v(:,2))