% defining the vdp eqn 
% y' = f(t), z'=f(t) is the system. Returns [y z]'
function result = vdp(t,r)
    mu = 100;
    result = zeros(2,1);
    result(1) = r(2) ;% y
    result(2) = mu*r(2)*(1-r(1)^2) - r(1) ; % z

