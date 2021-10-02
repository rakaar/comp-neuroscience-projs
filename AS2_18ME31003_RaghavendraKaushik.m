[x,y] = meshgrid(linspace(-2,2, 10), linspace(-2,2,10));

figure(1)
    hold on
        
        quiver(x,y, x-x^3 + y, x-y, 2, 'color',[0 0 0]);
        x_plot = linspace(-2,2,10)
        plot(x_plot,(x_plot.^3 -x_plot))
        plot(x_plot, x_plot)
    hold off
grid

figure(2)
    [x,y] = meshgrid(linspace(-5,5, 30), linspace(-5,5,30));
    quiver(x,y, x-x^3 + y, x-y, 2, 'color',[0 0 0]);
grid
